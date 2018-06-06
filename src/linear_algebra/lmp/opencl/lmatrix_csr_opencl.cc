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

/// @author Nico Trost, Benedikt Galler, Dimitar Lukarski, Martin Wlotzka

#include "../lmatrix_csr_cpu.h"
#include "opencl/lmatrix_csr_opencl.h"
#include "../lmp_mem.h"
#include "../init_vec_mat.h"
#include "../lmp_log.h"

#include <iostream>
#include <cassert>
#include <cstdlib>
#include <vector>

#ifdef WITH_OPENCL

#    ifdef __APPLE__
#        include <cl.h>
#    else
#        include <CL/cl.h>
#    endif

#    include "opencl/mem_opencl.h"
#    include "opencl/opencl_global.h"
#    include "opencl/opencl_kernel_mapper.h"
#    include "opencl/opencl_utils.h"

#else

#    define ERROR LOG_ERROR("ERROR no OPENCL support"); exit(-1);

#endif

// class OPENCL_CSR_lMatrix

template <typename ValueType>
OPENCL_CSR_lMatrix<ValueType>::OPENCL_CSR_lMatrix ( )
{
    this->platform_name_ = "OPENCL";
    this->platform_id_ = OPENCL;
    manager_initialized = false;
}

#ifdef WITH_OPENCL

template <typename ValueType>
OPENCL_CSR_lMatrix<ValueType>::OPENCL_CSR_lMatrix ( opencl_manager *man )
{

    this->platform_name_ = "OPENCL";
    this->platform_id_ = OPENCL;
    this->implementation_name_ = "OpenCL Stream";
    this->implementation_id_ = OPEN_CL;
    this->my_manager = new opencl_manager;
    my_manager = man;
    this->manager_initialized = true;

}
#endif

template <typename ValueType>
OPENCL_CSR_lMatrix<ValueType>::~OPENCL_CSR_lMatrix ( )
{
    this->Clear ( );
}

#ifdef WITH_OPENCL

template <typename ValueType>
void OPENCL_CSR_lMatrix<ValueType>::set_opencl_manager ( opencl_manager *man )
{

    if ( !manager_initialized ) this->my_manager = new opencl_manager;
    my_manager = man;
    this->manager_initialized = true;

}

template <typename ValueType>
opencl_manager* OPENCL_CSR_lMatrix<ValueType>::get_opencl_manager ( )
{

    assert ( manager_initialized );
    return this->my_manager;

}
#endif

template <typename ValueType>
void OPENCL_CSR_lMatrix<ValueType>::Clear ( )
{
#ifdef WITH_OPENCL

    if ( this->get_nnz ( ) > 0 )
    {
        opencl_memfreedev ( this->val );
        opencl_memfreedev ( this->col );
        opencl_memfreedev ( this->row );
    }

    this->nnz_ = 0;
    this->num_row_ = 0;
    this->num_col_ = 0;
    this->name_ = "";

#else
    ERROR;
#endif
}

template <typename ValueType>
void OPENCL_CSR_lMatrix<ValueType>::Zeros ( void )
{
#ifdef WITH_OPENCL

    assert ( manager_initialized );

    if ( this->get_nnz ( ) > 0 )
        openclzeros ( this->val, this->get_nnz ( ), my_manager->get_kernel ( 2 ), my_manager->get_command_queue ( ), this->global_threads, this->local_threads );

#else
    ERROR;
#endif
}

#ifdef WITH_OPENCL

template <typename ValueType>
void OPENCL_CSR_lMatrix<ValueType>::set_thread_blocks ( cl_device_id my_device, const int size )
{
    // set number of local threads to units * 8 (8 SPs per block)
    int thr = 128;
    local_threads = thr;
    global_threads = ( size / thr ) * thr;
    if ( size % thr != 0 ) global_threads += thr;
}
#endif

template <typename ValueType>
void OPENCL_CSR_lMatrix<ValueType>::Init ( const int init_nnz, const int init_num_row, const int init_num_col, const std::string init_name )
{
#ifdef WITH_OPENCL

    assert ( manager_initialized );
    assert ( init_nnz >= 0 );
    assert ( init_num_row >= 0 );
    assert ( init_num_col >= 0 );

    // delete the old structure
    if ( this->get_nnz ( ) > 0 )
        this->Clear ( );

    // allocate
    opencl_memallocdev<ValueType>( this->val, init_nnz, my_manager->get_context ( ) );

    opencl_memallocdev<int>( this->col, init_nnz, my_manager->get_context ( ) );

    opencl_memallocdev<int>( this->row, init_num_row + 1, my_manager->get_context ( ) );

    this->name_ = init_name;
    this->nnz_ = init_nnz;
    this->num_row_ = init_num_row;
    this->num_col_ = init_num_col;
    this->set_thread_blocks ( my_manager->get_device ( ), init_num_row );

    this->Zeros ( );
#else
    ERROR;
#endif
}

template <typename ValueType>
lMatrix<ValueType> &OPENCL_CSR_lMatrix<ValueType>::operator= ( const lMatrix<ValueType> &mat2 )
{
    if ( this == &mat2 )
        return *this;

    this->CopyFrom ( mat2 );
    return *this;
}

template <typename ValueType>
void OPENCL_CSR_lMatrix<ValueType>::CopyStructureFrom ( const lMatrix<ValueType> &mat2 )
{
#ifdef WITH_OPENCL

    assert ( manager_initialized );

    if ( this != &mat2 )
    {

        this->Init ( mat2.get_nnz ( ), mat2.get_num_row ( ), mat2.get_num_col ( ), mat2.get_name ( ) );

        // OPENCL = OPENCL
        if ( const OPENCL_CSR_lMatrix<ValueType> *casted_mat =
             dynamic_cast < const OPENCL_CSR_lMatrix<ValueType>* > ( &mat2 ) )
        {

            assert ( this->get_nnz ( ) == mat2.get_nnz ( ) );
            assert ( this->get_num_row ( ) == mat2.get_num_row ( ) );
            assert ( this->get_num_col ( ) == mat2.get_num_col ( ) );

            if ( this->get_nnz ( ) > 0 )
            {
                opencl_memcpydev<int>( this->col, casted_mat->col, mat2.get_nnz ( ), my_manager->get_command_queue ( ) );
                opencl_memcpydev<int>( this->row, casted_mat->row, mat2.get_num_row ( ) + 1, my_manager->get_command_queue ( ) );
            }

        }
        else
        {

            // OPENCL = CPU
            if ( const CPU_CSR_lMatrix<ValueType> *casted_mat2 =
                 dynamic_cast < const CPU_CSR_lMatrix<ValueType>* > ( &mat2 ) )
            {

                assert ( this->get_nnz ( ) == mat2.get_nnz ( ) );
                assert ( this->get_num_row ( ) == mat2.get_num_row ( ) );
                assert ( this->get_num_col ( ) == mat2.get_num_col ( ) );

                // copy
                if ( this->get_nnz ( ) > 0 )
                {
                    opencl_memcpy2dev<int>( this->col, casted_mat2->matrix.col, mat2.get_nnz ( ), my_manager->get_command_queue ( ) );
                    opencl_memcpy2dev<int>( this->row, casted_mat2->matrix.row, mat2.get_num_row ( ) + 1, my_manager->get_command_queue ( ) );
                }

            }
            else
            {

                // unsupported type
                LOG_ERROR ( "ERROR OPENCL_CSR_lMatrix<ValueType>::CopyStructureFrom unsupported matrix type" );
                this->print ( );
                mat2.print ( );
                exit ( -1 );
            }
        }
    }
#else
    ERROR;
#endif
}

template <typename ValueType>
void OPENCL_CSR_lMatrix<ValueType>::CopyFrom ( const lMatrix<ValueType> &mat2 )
{
#ifdef WITH_OPENCL

    assert ( manager_initialized );

    if ( this != &mat2 )
    {
        // OPENCL = OPENCL
        if ( const OPENCL_CSR_lMatrix<ValueType> *casted_mat =
             dynamic_cast < const OPENCL_CSR_lMatrix<ValueType>* > ( &mat2 ) )
        {

            assert ( this->get_nnz ( ) == mat2.get_nnz ( ) );
            assert ( this->get_num_row ( ) == mat2.get_num_row ( ) );
            assert ( this->get_num_col ( ) == mat2.get_num_col ( ) );

            if ( this->get_nnz ( ) > 0 )
            {
                opencl_memcpydev<ValueType>( this->val, casted_mat->val, mat2.get_nnz ( ), my_manager->get_command_queue ( ) );
            }
        }
        else
        {
            // OPENCL = CPU
            if ( const CPU_CSR_lMatrix<ValueType> *casted_mat2 =
                 dynamic_cast < const CPU_CSR_lMatrix<ValueType>* > ( &mat2 ) )
            {

                assert ( this->get_nnz ( ) == mat2.get_nnz ( ) );
                assert ( this->get_num_row ( ) == mat2.get_num_row ( ) );
                assert ( this->get_num_col ( ) == mat2.get_num_col ( ) );

                // copy
                if ( this->get_nnz ( ) > 0 )
                {
                    opencl_memcpy2dev<ValueType>( this->val, casted_mat2->matrix.val, mat2.get_nnz ( ), my_manager->get_command_queue ( ) );
                }
            }
            else
            {
                // unsupported type
                LOG_ERROR ( "ERROR OPENCL_CSR_lMatrix<ValueType>::operator= ; unsupported matrix type" );
                this->print ( );
                mat2.print ( );
                exit ( -1 );
            }
        }
    }
#else
    ERROR;
#endif
}

template <typename ValueType>
void OPENCL_CSR_lMatrix<ValueType>::CopyTo ( lMatrix<ValueType> &mat2 ) const
{
    if ( this != &mat2 )
    {

        // OPENCL CSR = OPENCL CSR
        if ( const OPENCL_CSR_lMatrix<ValueType> *casted_mat =
             dynamic_cast < const OPENCL_CSR_lMatrix<ValueType>* > ( &mat2 ) )
        {

            assert ( this->get_nnz ( ) == mat2.get_nnz ( ) );
            assert ( this->get_num_row ( ) == mat2.get_num_row ( ) );
            assert ( this->get_num_col ( ) == mat2.get_num_col ( ) );

            if ( this->get_nnz ( ) > 0 )
            {
#ifdef WITH_OPENCL
                opencl_memcpydev<ValueType>( casted_mat->val, this->val, this->get_nnz ( ), my_manager->get_command_queue ( ) );
#else
                ERROR;
#endif
            }

        }
        else
        {

            // OPENCL = CPU
            if ( const CPU_CSR_lMatrix<ValueType> *casted_mat2 =
                 dynamic_cast < const CPU_CSR_lMatrix<ValueType>* > ( &mat2 ) )
            {

                assert ( this->get_nnz ( ) == mat2.get_nnz ( ) );
                assert ( this->get_num_row ( ) == mat2.get_num_row ( ) );
                assert ( this->get_num_col ( ) == mat2.get_num_col ( ) );

                // copy
                if ( this->get_nnz ( ) > 0 )
                {
#ifdef WITH_OPENCL
                    opencl_memcpy2host<ValueType>( casted_mat2->matrix.val, this->val, this->get_nnz ( ), my_manager->get_command_queue ( ) );
#else
                    ERROR;
#endif
                }

            }
            else
            {

                LOG_ERROR ( "CPU_CSR_lMatrix<ValueType>::CopyTo unsupported matrix type" );
                this->print ( );
                mat2.print ( );
                exit ( -1 );

            }
        }
    }
}

template <typename ValueType>
void OPENCL_CSR_lMatrix<ValueType>::CopyStructureTo ( lMatrix<ValueType> &mat2 ) const
{
#ifdef WITH_OPENCL
    assert ( manager_initialized );
    if ( this != &mat2 )
    {

        mat2.Init ( this->get_nnz ( ), this->get_num_row ( ), this->get_num_col ( ), this->get_name ( ) );

        // OPENCL = OPENCL
        if ( const OPENCL_CSR_lMatrix<ValueType> *casted_mat =
             dynamic_cast < const OPENCL_CSR_lMatrix<ValueType>* > ( &mat2 ) )
        {

            assert ( this->get_nnz ( ) == mat2.get_nnz ( ) );
            assert ( this->get_num_row ( ) == mat2.get_num_row ( ) );
            assert ( this->get_num_col ( ) == mat2.get_num_col ( ) );

            if ( this->get_nnz ( ) > 0 )
            {
                opencl_memcpydev<int>( casted_mat->col, this->col, this->get_nnz ( ), my_manager->get_command_queue ( ) );
                opencl_memcpydev<int>( casted_mat->row, this->row, this->get_num_row ( ) + 1, my_manager->get_command_queue ( ) );
            }

        }
        else
        {

            // OPENCL = CPU
            if ( const CPU_CSR_lMatrix<ValueType> *casted_mat2 =
                 dynamic_cast < const CPU_CSR_lMatrix<ValueType>* > ( &mat2 ) )
            {

                assert ( this->get_nnz ( ) == mat2.get_nnz ( ) );
                assert ( this->get_num_row ( ) == mat2.get_num_row ( ) );
                assert ( this->get_num_col ( ) == mat2.get_num_col ( ) );

                // copy
                if ( this->get_nnz ( ) > 0 )
                {
                    opencl_memcpy2host<int>( casted_mat2->matrix.col, this->col, this->get_nnz ( ), my_manager->get_command_queue ( ) );
                    opencl_memcpy2host<int>( casted_mat2->matrix.row, this->row, this->get_num_row ( ) + 1, my_manager->get_command_queue ( ) );
                }

            }
            else
            {

                // unsupported type
                LOG_ERROR ( "ERROR OPENCL_CSR_lMatrix<ValueType>::CopyStructureFrom unsupported matrix type" );
                this->print ( );
                mat2.print ( );
                exit ( -1 );
            }
        }
    }
#else
    ERROR;
#endif
}

template <typename ValueType>
void OPENCL_CSR_lMatrix<ValueType>::ConvertFrom ( const lMatrix<ValueType> &mat2 )
{
    // unsupported type
    LOG_ERROR ( "ERROR OPENCL_CSR_lMatrix<ValueType>::ConvertFrom not implemented yet" );
    this->print ( );
    mat2.print ( );
    exit ( -1 );
}

template <typename ValueType>
lMatrix<ValueType> *OPENCL_CSR_lMatrix<ValueType>::extract_submatrix ( const int start_row, const int start_col,
                                                                       const int end_row, const int end_col ) const
{
#ifdef WITH_OPENCL
    lMatrix<ValueType> *cpu_matrix = new CPUsimple_CSR_lMatrix<ValueType>;
    lMatrix<ValueType> *cpu_sub_matrix;
    lMatrix<ValueType> *opencl_sub_matrix;

    cpu_matrix->CloneFrom ( *this );

    cpu_sub_matrix = cpu_matrix->extract_submatrix ( start_row, start_col,
                                                     end_row, end_col );

    std::string sub_mat_name;
    sub_mat_name = "sub matrix from ";
    sub_mat_name.append ( this->get_name ( ) );

    opencl_sub_matrix = new OPENCL_CSR_lMatrix<ValueType>( this->my_manager );
    opencl_sub_matrix->Init ( cpu_sub_matrix->get_nnz ( ), cpu_sub_matrix->get_num_row ( ), cpu_sub_matrix->get_num_col ( ), sub_mat_name );

    opencl_sub_matrix->CloneFrom ( *cpu_sub_matrix );

    delete cpu_matrix;
    delete cpu_sub_matrix;

    return opencl_sub_matrix;
#else
    ERROR;
#endif
}

template <typename ValueType>
void OPENCL_CSR_lMatrix<ValueType>::VectorMultAdd ( const OPENCL_lVector<ValueType> &invec, OPENCL_lVector<ValueType> *outvec ) const
{
#ifdef WITH_OPENCL

    assert ( invec .get_size ( ) >= 0 );
    assert ( outvec->get_size ( ) >= 0 );
    assert ( invec .get_size ( ) == this->get_num_col ( ) );
    assert ( ( outvec->get_size ( ) == this->get_num_row ( ) ) ||
             ( invec .get_size ( ) == 0 ) );
    assert ( manager_initialized );

    if ( this->get_nnz ( ) > 0 )
        openclcsrvectormult ( this->val, this->col, this->row, invec.buffer, outvec->buffer, this->get_num_row ( ), my_manager->get_kernel ( 9 ), my_manager->get_command_queue ( ), this->global_threads, this->local_threads );

#else
    ERROR;
#endif
}

template <typename ValueType>
void OPENCL_CSR_lMatrix<ValueType>::VectorMultAdd ( const lVector<ValueType> &invec, lVector<ValueType> *outvec ) const
{
#ifdef WITH_OPENCL

    assert ( invec .get_size ( ) >= 0 );
    assert ( outvec->get_size ( ) >= 0 );
    assert ( invec .get_size ( ) == this->get_num_col ( ) );
    assert ( ( outvec->get_size ( ) == this->get_num_row ( ) ) ||
             ( invec .get_size ( ) == 0 ) );
    assert ( manager_initialized );

    const OPENCL_lVector<ValueType> *casted_invec = dynamic_cast < const OPENCL_lVector<ValueType>* > ( &invec );
    OPENCL_lVector<ValueType> *casted_outvec = dynamic_cast < OPENCL_lVector<ValueType>* > ( outvec );

    if ( ( casted_invec == NULL ) && ( casted_outvec == NULL ) )
    {
        LOG_ERROR ( "ERROR OPENCL_CSR_lMatrix<ValueType>::VectorMultAdd unsupported in or out vector" );
        this->print ( );
        invec.print ( );
        outvec->print ( );
        exit ( -1 );
    }

    this->VectorMultAdd ( *casted_invec, casted_outvec );
#else
    ERROR;
#endif
}

template <typename ValueType>
void OPENCL_CSR_lMatrix<ValueType>::VectorMult ( const OPENCL_lVector<ValueType> &invec, OPENCL_lVector<ValueType> *outvec ) const
{
#ifdef WITH_OPENCL

    assert ( invec .get_size ( ) >= 0 );
    assert ( outvec->get_size ( ) >= 0 );
    assert ( invec .get_size ( ) == this->get_num_col ( ) );
    assert ( outvec->get_size ( ) == this->get_num_row ( ) );
    assert ( manager_initialized );

    if ( this->get_nnz ( ) > 0 )
        openclcsrvectormult ( this->val, this->col, this->row, invec.buffer, outvec->buffer, this->get_num_row ( ), my_manager->get_kernel ( 8 ), my_manager->get_command_queue ( ), this->global_threads, this->local_threads );

#else
    ERROR;
#endif
}

template <typename ValueType>
void OPENCL_CSR_lMatrix<ValueType>::VectorMult ( const lVector<ValueType> &invec, lVector<ValueType> *outvec ) const
{
#ifdef WITH_OPENCL

    assert ( this->get_num_col ( ) == invec.get_size ( ) );
    assert ( this->get_num_row ( ) == outvec->get_size ( ) );
    assert ( manager_initialized );

    const OPENCL_lVector<ValueType> *casted_invec = dynamic_cast < const OPENCL_lVector<ValueType>* > ( &invec );
    OPENCL_lVector<ValueType> *casted_outvec = dynamic_cast < OPENCL_lVector<ValueType>* > ( outvec );

    if ( ( casted_invec == NULL ) && ( casted_outvec == NULL ) )
    {
        LOG_ERROR ( "ERROR OPENCL_CSR_lMatrix<ValueType>::VectorMult unsupported in or out vector" );
        this->print ( );
        invec.print ( );
        outvec->print ( );
        exit ( -1 );
    }

    this->VectorMult ( *casted_invec, casted_outvec );

#else
    ERROR;
#endif
}

template <typename ValueType>
lMatrix<ValueType> *OPENCL_CSR_lMatrix<ValueType>::CloneWithoutContent ( ) const
{
#ifdef WITH_OPENCL

    std::string cloned_name;
    cloned_name = "clone from ";
    cloned_name.append ( this->name_ );

    lMatrix<ValueType> *new_matrix = new OPENCL_CSR_lMatrix<ValueType>( this->my_manager );
    new_matrix->Init ( this->get_nnz ( ), this->get_num_row ( ), this->get_num_col ( ), cloned_name );
    new_matrix->Zeros ( );
    return new_matrix;

#else
    ERROR;
#endif
}

template <typename ValueType>
void OPENCL_CSR_lMatrix<ValueType>::CastFrom ( const lMatrix<double>& other )
{
    assert ( this->get_num_row ( ) == other.get_num_row ( ) );
    assert ( this->get_num_col ( ) == other.get_num_col ( ) );
    assert ( this->get_nnz ( ) == other.get_nnz ( ) );

    LOG_ERROR ( "OPENCL_CSR_lMatrix<ValueType>::CastFrom<double> not yet implemented." );
    this->print ( );
    other.print ( );
    exit ( -1 );
}

template <typename ValueType>
void OPENCL_CSR_lMatrix<ValueType>::CastFrom ( const lMatrix<float>& other )
{
    assert ( this->get_num_row ( ) == other.get_num_row ( ) );
    assert ( this->get_num_col ( ) == other.get_num_col ( ) );
    assert ( this->get_nnz ( ) == other.get_nnz ( ) );

    LOG_ERROR ( "OPENCL_CSR_lMatrix<ValueType>::CastFrom<float> not yet implemented." );
    this->print ( );
    other.print ( );
    exit ( -1 );
}

template <typename ValueType>
void OPENCL_CSR_lMatrix<ValueType>::CastTo ( lMatrix<double>& other ) const
{
    assert ( this->get_num_row ( ) == other.get_num_row ( ) );
    assert ( this->get_num_col ( ) == other.get_num_col ( ) );
    assert ( this->get_nnz ( ) == other.get_nnz ( ) );

    LOG_ERROR ( "OPENCL_CSR_lMatrix<ValueType>::CastTo<double> not yet implemented." );
    this->print ( );
    other.print ( );
    exit ( -1 );
}

template <typename ValueType>
void OPENCL_CSR_lMatrix<ValueType>::CastTo ( lMatrix<float>& other ) const
{
    assert ( this->get_num_row ( ) == other.get_num_row ( ) );
    assert ( this->get_num_col ( ) == other.get_num_col ( ) );
    assert ( this->get_nnz ( ) == other.get_nnz ( ) );

    LOG_ERROR ( "OPENCL_CSR_lMatrix<ValueType>::CastTo<float> not yet implemented." );
    this->print ( );
    other.print ( );
    exit ( -1 );
}

template <typename ValueType>
void OPENCL_CSR_lMatrix<ValueType>::SwapDiagElementsToRowFront ( void )
{
    LOG_ERROR ( "OPENCL_CSR_lMatrix<ValueType>::SwapDiagElementsToRowFront() not yet implemented." );
    this->print ( );
    exit ( -1 );
}

template <typename ValueType>
void OPENCL_CSR_lMatrix<ValueType>::VectorMultNoDiag ( const lVector<ValueType> &in,
                                                       lVector<ValueType> *out ) const
{
    LOG_ERROR ( "OPENCL_CSR_lMatrix<ValueType>::VectorMultNoDiag not yet implemented." );
    this->print ( );
    exit ( -1 );
}

template class OPENCL_CSR_lMatrix<double>;
template class OPENCL_CSR_lMatrix<float>;
