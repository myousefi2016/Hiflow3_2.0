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

#include "cuda/lmatrix_coo_gpu.h"
#include "../lmp_mem.h"
#include "../lmp_log.h"

#include <iostream>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

using namespace hiflow::la;

// class GPU_COO_lMatrix

template <typename ValueType>
GPU_COO_lMatrix<ValueType>::GPU_COO_lMatrix ( )
{
    this->platform_name_ = "GPU (CUDA)";
    this->platform_id_ = GPU;

    this->matrix.val = NULL;
    this->matrix.col = NULL;
    this->matrix.row = NULL;

    this->thread_block_size_ = 0;
    this->thread_block_ = 0;
}

template <typename ValueType>
GPU_COO_lMatrix<ValueType>::~GPU_COO_lMatrix ( )
{
    this->Clear ( );
}

template <typename ValueType>
void GPU_COO_lMatrix<ValueType>::Clear ( )
{
    if ( this->matrix.val != NULL )
        memfreedev ( this->matrix.val );
    this->matrix.val = NULL;

    if ( this->matrix.col != NULL )
        memfreedev ( this->matrix.col );
    this->matrix.col = NULL;

    if ( this->matrix.row != NULL )
        memfreedev ( this->matrix.row );
    this->matrix.row = NULL;

    this->nnz_ = 0;
    this->num_row_ = 0;
    this->num_col_ = 0;
    this->name_ = "";
}

template <typename ValueType>
void GPU_COO_lMatrix<ValueType>::Zeros ( void )
{
    memsetdev ( this->matrix.val, 0, this->get_nnz ( ), sizeof (ValueType ) );
}

template <typename ValueType>
void GPU_COO_lMatrix<ValueType>::set_thread_blocks ( int thread_block, int thread_block_size )
{
    this->thread_block_ = thread_block;
    this->thread_block_size_ = thread_block_size;
}

template <typename ValueType>
void GPU_COO_lMatrix<ValueType>::set_thread_blocks ( void )
{
    // default values
    int thread_block_size = 256;
    int thread_block = this->get_num_row ( ) / thread_block_size;

    if ( thread_block_size * thread_block < this->get_num_row ( ) )
        thread_block++;

    this->set_thread_blocks ( thread_block, thread_block_size );
}

template <typename ValueType>
void GPU_COO_lMatrix<ValueType>::Init ( const int init_nnz, const int init_num_row, const int init_num_col, const std::string init_name )
{

    assert ( init_nnz >= 0 );
    assert ( init_num_row >= 0 );
    assert ( init_num_col >= 0 );

    // delete the old structure
    this->Clear ( );

    // allocate
    memallocdev ( ( void** ) & this->matrix.val, init_nnz, sizeof (ValueType ) );
    memsetdev ( this->matrix.val, 0, init_nnz, sizeof (ValueType ) );

    memallocdev ( ( void** ) & this->matrix.col, init_nnz, sizeof (int ) );
    memsetdev ( this->matrix.col, 0, init_nnz, sizeof (int ) );

    memallocdev ( ( void** ) & this->matrix.row, init_nnz, sizeof (int ) );
    memsetdev ( this->matrix.row, 0, init_nnz, sizeof (int ) );

    this->name_ = init_name;
    this->nnz_ = init_nnz;
    this->num_row_ = init_num_row;
    this->num_col_ = init_num_col;
    this->set_thread_blocks ( );

    this->Zeros ( );

}

template <typename ValueType>
lMatrix<ValueType> &GPU_COO_lMatrix<ValueType>::operator= ( const lMatrix<ValueType> &mat2 )
{

    if ( this == &mat2 )
        return *this;

    this->CopyFrom ( mat2 );
    return *this;
}

template <typename ValueType>
void GPU_COO_lMatrix<ValueType>::CopyFrom ( const lMatrix<ValueType> &mat2 )
{

    if ( this != &mat2 )
    {

        // GPU = GPU
        if ( const GPU_COO_lMatrix<ValueType> *casted_mat =
             dynamic_cast < const GPU_COO_lMatrix<ValueType>* > ( &mat2 ) )
        {

            assert ( this->get_nnz ( ) == mat2.get_nnz ( ) );
            assert ( this->get_num_row ( ) == mat2.get_num_row ( ) );
            assert ( this->get_num_col ( ) == mat2.get_num_col ( ) );

            // copy
            if ( this->get_nnz ( ) > 0 )
            {
                memcpydev ( this->matrix.val, casted_mat->matrix.val, this->get_nnz ( ) );
            }

        }
        else
        {

            // GPU = CPU
            if ( const CPU_COO_lMatrix<ValueType> *casted_mat =
                 dynamic_cast < const CPU_COO_lMatrix<ValueType>* > ( &mat2 ) )
            {

                assert ( this->get_nnz ( ) == mat2.get_nnz ( ) );
                assert ( this->get_num_row ( ) == mat2.get_num_row ( ) );
                assert ( this->get_num_col ( ) == mat2.get_num_col ( ) );

                // copy
                if ( this->get_nnz ( ) > 0 )
                {
                    memcpy2dev ( this->matrix.val, casted_mat->matrix.val, this->get_nnz ( ) );
                }

            }
            else
            {

                // unsupported type
                LOG_ERROR ( "ERROR GPU_COO_lMatrix<ValueType>::CopyFrom() unsupported matrix type" );
                this->print ( );
                mat2.print ( );
                exit ( -1 );
            }
        }
    }
}

template <typename ValueType>
void GPU_COO_lMatrix<ValueType>::CopyTo ( lMatrix<ValueType> &mat2 ) const
{
    if ( this != &mat2 )
    {

        // GPU COO = GPU COO
        if ( const GPU_COO_lMatrix<ValueType> *casted_mat =
             dynamic_cast < const GPU_COO_lMatrix<ValueType>* > ( &mat2 ) )
        {

            assert ( this->get_nnz ( ) == mat2.get_nnz ( ) );
            assert ( this->get_num_row ( ) == mat2.get_num_row ( ) );
            assert ( this->get_num_col ( ) == mat2.get_num_col ( ) );

            if ( this->get_nnz ( ) > 0 )
            {
                memcpydev ( casted_mat->matrix.val, this->matrix.val, this->get_nnz ( ) );
            }

        }
        else
        {

            // GPU = CPU
            if ( const CPU_COO_lMatrix<ValueType> *casted_mat =
                 dynamic_cast < const CPU_COO_lMatrix<ValueType>* > ( &mat2 ) )
            {

                assert ( this->get_nnz ( ) == mat2.get_nnz ( ) );
                assert ( this->get_num_row ( ) == mat2.get_num_row ( ) );
                assert ( this->get_num_col ( ) == mat2.get_num_col ( ) );

                // copy
                if ( this->get_nnz ( ) > 0 )
                {
                    memcpy2host ( casted_mat->matrix.val, this->matrix.val, this->get_nnz ( ) );
                }

            }
            else
            {

                LOG_ERROR ( "CPU_COO_lMatrix<ValueType>::CopyTo unsupported matrix type" );
                this->print ( );
                mat2.print ( );
                exit ( -1 );

            }
        }
    }
}

template <typename ValueType>
void GPU_COO_lMatrix<ValueType>::CopyStructureFrom ( const lMatrix<ValueType> &mat2 )
{

    if ( this != &mat2 )
    {

        this->Init ( mat2.get_nnz ( ), mat2.get_num_row ( ), mat2.get_num_col ( ), mat2.get_name ( ) );

        // GPU COO = GPU COO
        if ( const GPU_COO_lMatrix<ValueType> *casted_mat =
             dynamic_cast < const GPU_COO_lMatrix<ValueType>* > ( &mat2 ) )
        {

            assert ( this->get_nnz ( ) == mat2.get_nnz ( ) );
            assert ( this->get_num_row ( ) == mat2.get_num_row ( ) );
            assert ( this->get_num_col ( ) == mat2.get_num_col ( ) );

            // copy
            if ( this->get_nnz ( ) > 0 )
            {
                memcpydev ( this->matrix.col, casted_mat->matrix.col, this->get_nnz ( ) );
                memcpydev ( this->matrix.row, casted_mat->matrix.row, this->get_nnz ( ) );
            }

        }
        else
        {

            // GPU COO = CPU COO
            if ( const CPU_COO_lMatrix<ValueType> *casted_mat =
                 dynamic_cast < const CPU_COO_lMatrix<ValueType>* > ( &mat2 ) )
            {

                assert ( this->get_nnz ( ) == mat2.get_nnz ( ) );
                assert ( this->get_num_row ( ) == mat2.get_num_row ( ) );
                assert ( this->get_num_col ( ) == mat2.get_num_col ( ) );

                // copy
                if ( this->get_nnz ( ) > 0 )
                {
                    memcpy2host ( this->matrix.col, casted_mat->matrix.col, this->get_nnz ( ) );
                    memcpy2host ( this->matrix.row, casted_mat->matrix.row, this->get_nnz ( ) );
                }

            }
            else
            {

                // unsupported type
                LOG_ERROR ( "GPU_COO_lMatrix<ValueType>::CopyStructureFrom unsupported matrix type" );
                this->print ( );
                mat2.print ( );
                exit ( -1 );
            }
        }
    }
}

template <typename ValueType>
void GPU_COO_lMatrix<ValueType>::CopyStructureTo ( lMatrix<ValueType> &mat2 ) const
{

    if ( this != &mat2 )
    {

        mat2.Init ( this->get_nnz ( ), this->get_num_row ( ), this->get_num_col ( ), this->get_name ( ) );

        // GPU COO = GPU COO
        if ( const GPU_COO_lMatrix<ValueType> *casted_mat =
             dynamic_cast < const GPU_COO_lMatrix<ValueType>* > ( &mat2 ) )
        {

            assert ( this->get_nnz ( ) == mat2.get_nnz ( ) );
            assert ( this->get_num_row ( ) == mat2.get_num_row ( ) );
            assert ( this->get_num_col ( ) == mat2.get_num_col ( ) );

            // copy
            if ( this->get_nnz ( ) > 0 )
            {
                memcpydev ( casted_mat->matrix.col, this->matrix.col, this->get_nnz ( ) );
                memcpydev ( casted_mat->matrix.row, this->matrix.row, this->get_nnz ( ) );
            }

        }
        else
        {

            // CPU COO = GPU COO
            if ( const CPU_COO_lMatrix<ValueType> *casted_mat =
                 dynamic_cast < const CPU_COO_lMatrix<ValueType>* > ( &mat2 ) )
            {

                assert ( this->get_nnz ( ) == mat2.get_nnz ( ) );
                assert ( this->get_num_row ( ) == mat2.get_num_row ( ) );
                assert ( this->get_num_col ( ) == mat2.get_num_col ( ) );

                if ( this->get_nnz ( ) > 0 )
                {
                    memcpy2host ( casted_mat->matrix.col, this->matrix.col, this->get_nnz ( ) );
                    memcpy2host ( casted_mat->matrix.row, this->matrix.row, this->get_nnz ( ) );
                }

            }
            else
            {

                LOG_ERROR ( "GPU_COO_lMatrix<ValueType>::CopyStructureTo unsupported matrix type" );
                this->print ( );
                mat2.print ( );
                exit ( -1 );

            }
        }
    }
}

template <typename ValueType>
void GPU_COO_lMatrix<ValueType>::ConvertFrom ( const lMatrix<ValueType> &mat2 )
{
    // unsupported type
    LOG_ERROR ( "ERROR GPU_COO_lMatrix<ValueType>::ConvertFrom unsupported matrix type" );
    this->print ( );
    mat2.print ( );
    exit ( -1 );
}

template <typename ValueType>
void GPU_COO_lMatrix<ValueType>::Sync ( void ) const
{
    cuda_sync_threads ( );
}

template <typename ValueType>
void GPU_COO_lMatrix<ValueType>::CastFrom ( const hiflow::la::lMatrix<double>& other )
{
    assert ( this->get_num_row ( ) == other.get_num_row ( ) );
    assert ( this->get_num_col ( ) == other.get_num_col ( ) );
    assert ( this->get_nnz ( ) == other.get_nnz ( ) );

    if ( const COO_lMatrix<double>* other_coo
         = dynamic_cast < const COO_lMatrix<double>* > ( &other ) )
    {
        // other matrix is also in COO format
        if ( const GPU_COO_lMatrix<double>* other_gpu_coo
             = dynamic_cast < const GPU_COO_lMatrix<double>* > ( other_coo ) )
        {
            // GPU from GPU
            this->CastFrom ( *other_gpu_coo );
        }
        else if ( const CPU_COO_lMatrix<double>* other_cpu_coo
                  = dynamic_cast < const CPU_COO_lMatrix<double>* > ( other_coo ) )
        {
            // GPU from CPU
            this->CastFrom ( *other_cpu_coo );
        }
        else
        {
            LOG_ERROR ( "GPU_COO_lMatrix<ValueType>::CastFrom<double> called with non-GPU and non-CPU matrix argument." );
            this->print ( );
            other.print ( );
            exit ( -1 );
        }
    }
    else
    {
        LOG_ERROR ( "GPU_COO_lMatrix<ValueType>::CastFrom<double> called with non-COO matrix argument." );
        this->print ( );
        other.print ( );
        exit ( -1 );
    }
}

template <typename ValueType>
void GPU_COO_lMatrix<ValueType>::CastFrom ( const hiflow::la::lMatrix<float>& other )
{
    assert ( this->get_num_row ( ) == other.get_num_row ( ) );
    assert ( this->get_num_col ( ) == other.get_num_col ( ) );
    assert ( this->get_nnz ( ) == other.get_nnz ( ) );

    if ( const COO_lMatrix<float>* other_coo
         = dynamic_cast < const COO_lMatrix<float>* > ( &other ) )
    {
        // other matrix is also in COO format
        if ( const GPU_COO_lMatrix<float>* other_gpu_coo
             = dynamic_cast < const GPU_COO_lMatrix<float>* > ( other_coo ) )
        {
            // GPU from GPU
            this->CastFrom ( *other_gpu_coo );
        }
        else if ( const CPU_COO_lMatrix<float>* other_cpu_coo
                  = dynamic_cast < const CPU_COO_lMatrix<float>* > ( other_coo ) )
        {
            // GPU from CPU
            this->CastFrom ( *other_cpu_coo );
        }
        else
        {
            LOG_ERROR ( "GPU_COO_lMatrix<ValueType>::CastFrom<float> called with non-GPU and non-CPU matrix argument." );
            this->print ( );
            other.print ( );
            exit ( -1 );
        }
    }
    else
    {
        LOG_ERROR ( "GPU_COO_lMatrix<ValueType>::CastFrom<float> called with non-COO matrix argument." );
        this->print ( );
        other.print ( );
        exit ( -1 );
    }
}

template <typename ValueType>
void GPU_COO_lMatrix<ValueType>::CastFrom ( const GPU_COO_lMatrix<double>& other )
{
    assert ( this->get_num_row ( ) == other.get_num_row ( ) );
    assert ( this->get_num_col ( ) == other.get_num_col ( ) );
    assert ( this->get_nnz ( ) == other.get_nnz ( ) );

    cudacastfromdouble ( this->get_nnz ( ), this->matrix.val, other.matrix.val, this->thread_block_size_ );
}

template <typename ValueType>
void GPU_COO_lMatrix<ValueType>::CastFrom ( const GPU_COO_lMatrix<float>& other )
{
    assert ( this->get_num_row ( ) == other.get_num_row ( ) );
    assert ( this->get_num_col ( ) == other.get_num_col ( ) );
    assert ( this->get_nnz ( ) == other.get_nnz ( ) );

    cudacastfromfloat ( this->get_nnz ( ), this->matrix.val, other.matrix.val, this->thread_block_size_ );
}

template <typename ValueType>
void GPU_COO_lMatrix<ValueType>::CastFrom ( const CPU_COO_lMatrix<double>& other )
{
    assert ( this->get_num_row ( ) == other.get_num_row ( ) );
    assert ( this->get_num_col ( ) == other.get_num_col ( ) );
    assert ( this->get_nnz ( ) == other.get_nnz ( ) );

    ValueType *vals = new ValueType[this->get_nnz ( )];

    for ( int i = 0; i < this->get_nnz ( ); ++i )
    {
        vals[i] = static_cast < ValueType > ( other.matrix.val[i] );
    }

    memcpy2dev ( this->matrix.val, vals, this->get_nnz ( ) );

    delete [] vals;
}

template <typename ValueType>
void GPU_COO_lMatrix<ValueType>::CastFrom ( const CPU_COO_lMatrix<float>& other )
{
    assert ( this->get_num_row ( ) == other.get_num_row ( ) );
    assert ( this->get_num_col ( ) == other.get_num_col ( ) );
    assert ( this->get_nnz ( ) == other.get_nnz ( ) );

    ValueType *vals = new ValueType[this->get_nnz ( )];

    for ( int i = 0; i < this->get_nnz ( ); ++i )
    {
        vals[i] = static_cast < ValueType > ( other.matrix.val[i] );
    }

    memcpy2dev ( this->matrix.val, vals, this->get_nnz ( ) );

    delete [] vals;
}

template <typename ValueType>
void GPU_COO_lMatrix<ValueType>::CastTo ( hiflow::la::lMatrix<double>& other ) const
{
    assert ( this->get_num_row ( ) == other.get_num_row ( ) );
    assert ( this->get_num_col ( ) == other.get_num_col ( ) );
    assert ( this->get_nnz ( ) == other.get_nnz ( ) );

    if ( COO_lMatrix<double>* other_coo
         = dynamic_cast < COO_lMatrix<double>* > ( &other ) )
    {
        // other matrix is also in COO format
        if ( GPU_COO_lMatrix<double>* other_gpu_coo
             = dynamic_cast < GPU_COO_lMatrix<double>* > ( other_coo ) )
        {
            // GPU to GPU
            this->CastTo ( *other_gpu_coo );
        }
        else if ( CPU_COO_lMatrix<double>* other_cpu_coo
                  = dynamic_cast < CPU_COO_lMatrix<double>* > ( other_coo ) )
        {
            // GPU to CPU
            this->CastTo ( *other_cpu_coo );
        }
        else
        {
            LOG_ERROR ( "GPU_COO_lMatrix<ValueType>::CastTo<double> called with non-GPU and non-CPU matrix argument." );
            this->print ( );
            other.print ( );
            exit ( -1 );
        }
    }
    else
    {
        LOG_ERROR ( "GPU_COO_lMatrix<ValueType>::CastTo<double> called with non-COO matrix argument." );
        this->print ( );
        other.print ( );
        exit ( -1 );
    }
}

template <typename ValueType>
void GPU_COO_lMatrix<ValueType>::CastTo ( hiflow::la::lMatrix<float>& other ) const
{
    assert ( this->get_num_row ( ) == other.get_num_row ( ) );
    assert ( this->get_num_col ( ) == other.get_num_col ( ) );
    assert ( this->get_nnz ( ) == other.get_nnz ( ) );

    if ( COO_lMatrix<float>* other_coo
         = dynamic_cast < COO_lMatrix<float>* > ( &other ) )
    {
        // other matrix is also in COO format
        if ( GPU_COO_lMatrix<float>* other_gpu_coo
             = dynamic_cast < GPU_COO_lMatrix<float>* > ( other_coo ) )
        {
            // GPU to GPU
            this->CastTo ( *other_gpu_coo );
        }
        else if ( CPU_COO_lMatrix<float>* other_cpu_coo
                  = dynamic_cast < CPU_COO_lMatrix<float>* > ( other_coo ) )
        {
            // GPU to CPU
            this->CastTo ( *other_cpu_coo );
        }
        else
        {
            LOG_ERROR ( "GPU_COO_lMatrix<ValueType>::CastTo<float> called with non-GPU and non-CPU matrix argument." );
            this->print ( );
            other.print ( );
            exit ( -1 );
        }
    }
    else
    {
        LOG_ERROR ( "GPU_COO_lMatrix<ValueType>::CastTo<float> called with non-COO matrix argument." );
        this->print ( );
        other.print ( );
        exit ( -1 );
    }
}

template <typename ValueType>
void GPU_COO_lMatrix<ValueType>::CastTo ( GPU_COO_lMatrix<double>& other ) const
{
    assert ( this->get_num_row ( ) == other.get_num_row ( ) );
    assert ( this->get_num_col ( ) == other.get_num_col ( ) );
    assert ( this->get_nnz ( ) == other.get_nnz ( ) );

    cudacasttodouble ( this->get_nnz ( ), other.matrix.val, this->matrix.val, this->thread_block_size_ );
}

template <typename ValueType>
void GPU_COO_lMatrix<ValueType>::CastTo ( GPU_COO_lMatrix<float>& other ) const
{
    assert ( this->get_num_row ( ) == other.get_num_row ( ) );
    assert ( this->get_num_col ( ) == other.get_num_col ( ) );
    assert ( this->get_nnz ( ) == other.get_nnz ( ) );

    cudacasttofloat ( this->get_nnz ( ), other.matrix.val, this->matrix.val, this->thread_block_size_ );
}

template <typename ValueType>
void GPU_COO_lMatrix<ValueType>::CastTo ( CPU_COO_lMatrix<double>& other ) const
{
    assert ( this->get_num_row ( ) == other.get_num_row ( ) );
    assert ( this->get_num_col ( ) == other.get_num_col ( ) );
    assert ( this->get_nnz ( ) == other.get_nnz ( ) );

    ValueType *vals = new ValueType[this->get_nnz ( )];

    memcpy2host ( vals, this->matrix.val, this->get_nnz ( ) );

    for ( int i = 0; i < this->get_nnz ( ); ++i )
    {
        other.matrix.val[i] = static_cast < double > ( vals[i] );
    }

    delete [] vals;
}

template <typename ValueType>
void GPU_COO_lMatrix<ValueType>::CastTo ( CPU_COO_lMatrix<float>& other ) const
{
    assert ( this->get_num_row ( ) == other.get_num_row ( ) );
    assert ( this->get_num_col ( ) == other.get_num_col ( ) );
    assert ( this->get_nnz ( ) == other.get_nnz ( ) );

    ValueType *vals = new ValueType[this->get_nnz ( )];

    memcpy2host ( vals, this->matrix.val, this->get_nnz ( ) );

    for ( int i = 0; i < this->get_nnz ( ); ++i )
    {
        other.matrix.val[i] = static_cast < float > ( vals[i] );
    }

    delete [] vals;
}

template class GPU_COO_lMatrix<float>;
template class GPU_COO_lMatrix<double>;
