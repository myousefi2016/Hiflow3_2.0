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

#include "lmatrix.h"
#include "lmatrix_csr_cpu.h"
#include "init_vec_mat.h"
#include "lmp_log.h"

#include <iostream>
#include <stdlib.h>

namespace hiflow
{
    namespace la
    {

        template <typename ValueType>
        lMatrix<ValueType>::lMatrix ( )
        {
            this->nnz_ = 0;
            this->num_col_ = 0;
            this->num_row_ = 0;

            this->name_ = "";
            this->platform_name_ = "";
            this->platform_id_ = CPU;
            this->implementation_name_ = "";
            this->implementation_id_ = NAIVE;
            this->matrix_format_name_ = "";
            this->matrix_format_id_ = CSR; // default format

            this->symmetric_ = false;
        }

        template <typename ValueType>
        lMatrix<ValueType>::~lMatrix ( )
        {
        }

        template <typename ValueType>
        void lMatrix<ValueType>::get_as_coo ( std::vector<ValueType>& vals,
                                              std::vector<int>& rows,
                                              std::vector<int>& cols ) const
        {
            LOG_ERROR ( "lMatrix::get_as_coo not implemented for this matrix type." );
            this->print ( );
            exit ( -1 );
        }

        template <typename ValueType>
        lMatrix<ValueType> *lMatrix<ValueType>::CloneWithoutContent ( ) const
        {
            std::string cloned_name;
            cloned_name = "clone from ";
            cloned_name.append ( this->name_ );

            lMatrix<ValueType> *new_matrix = init_matrix<ValueType>( this->get_nnz ( ),
                    this->get_num_row ( ),
                    this->get_num_col ( ),
                    cloned_name,
                    this->get_platform ( ),
                    this->get_implementation ( ),
                    this->get_matrix_format ( ) );
            new_matrix->Zeros ( );
            return new_matrix;
        }

        template<typename ValueType>
        void lMatrix<ValueType>::Compress ( )
        {
        }

        template <typename ValueType>
        void lMatrix<ValueType>::print ( std::ostream &out ) const
        {

            LOG_INFO ( "lMatrix",
                       "name='" << this->get_name ( ) << "' Matrix format=" << this->get_matrix_format_name ( )
                       << " nnz=" << this->get_nnz ( ) << ", num_row=" << this->get_num_row ( )
                       << ", num_col=" << this->get_num_col ( )
                       << ", Precision=" << sizeof (ValueType )*8 << "bit"
                       << ", Platform:" << this->get_platform_name ( )
                       << ", Implementation:" << this->get_implementation_name ( )
                       << ", Symmetric check:" << this->symmetric_
                       );
        }

        template <typename ValueType>
        enum MATRIX_FORMAT lMatrix<ValueType>::get_matrix_format ( void ) const
        {
            return this->matrix_format_id_;
        }

        template <typename ValueType>
        std::string lMatrix<ValueType>::get_matrix_format_name ( void ) const
        {
            return this->matrix_format_name_;
        }

        template <typename ValueType>
        std::string lMatrix<ValueType>::get_name ( void ) const
        {
            return this->name_;
        }

        template <typename ValueType>
        enum IMPLEMENTATION lMatrix<ValueType>::get_implementation ( void ) const
        {
            return this->implementation_id_;
        }

        template <typename ValueType>
        std::string lMatrix<ValueType>::get_implementation_name ( void ) const
        {
            return this->implementation_name_;
        }

        template <typename ValueType>
        enum PLATFORM lMatrix<ValueType>::get_platform ( void ) const
        {
            return this->platform_id_;
        }

        template <typename ValueType>
        std::string lMatrix<ValueType>::get_platform_name ( void ) const
        {
            return this->platform_name_;
        }

        template <typename ValueType>
        void lMatrix<ValueType>::Sync ( void ) const
        {
            // do nothing
        }

        template <typename ValueType>
        void lMatrix<ValueType>::ReadFile ( const char* filename )
        {
            LOG_ERROR ( "lMatrix::ReadFile() does not support this matrix format" );
            this->print ( );
            exit ( -1 );
        }

        template <typename ValueType>
        void lMatrix<ValueType>::WriteFile ( const char* filename ) const
        {
            LOG_ERROR ( "lMatrix::WriteFile() does not support this matrix format" );
            this->print ( );
            exit ( -1 );
        }

        template <typename ValueType>
        void lMatrix<ValueType>::Pjacobi ( const lVector<ValueType> &invec, lVector<ValueType> *outvec ) const
        {
            LOG_ERROR ( "Preconditioner Jacobi does not support this matrix" );
            this->print ( );
            exit ( -1 );
        }

        template <typename ValueType>
        void lMatrix<ValueType>::Pgauss_seidel ( const lVector<ValueType> &invec, lVector<ValueType> *outvec ) const
        {
            LOG_ERROR ( "Preconditioner Gauss-Seidel does not support this matrix" );
            this->print ( );
            exit ( -1 );
        }

        template <typename ValueType>
        void lMatrix<ValueType>::Psor ( const ValueType omega, const lVector<ValueType> &invec, lVector<ValueType> *outvec ) const
        {
            LOG_ERROR ( "Preconditioner SOR does not support this matrix" );
            this->print ( );
            exit ( -1 );
        }

        template <typename ValueType>
        void lMatrix<ValueType>::Pssor ( const ValueType omega, const lVector<ValueType> &invec, lVector<ValueType> *outvec ) const
        {
            LOG_ERROR ( "Preconditioner SSOR does not support this matrix" );
            this->print ( );
            exit ( -1 );
        }

        template <typename ValueType>
        void lMatrix<ValueType>::Psgauss_seidel ( const lVector<ValueType> &invec, lVector<ValueType> *outvec ) const
        {
            LOG_ERROR ( "Preconditioner Symmetric Gauss Seidel does not support this matrix" );
            this->print ( );
            exit ( -1 );
        }

        template <typename ValueType>
        void lMatrix<ValueType>::BlockPsgauss_seidel ( const lVector<ValueType> &invec, lVector<ValueType> *outvec,
                                                       const int start_i, const int end_i ) const
        {
            LOG_ERROR ( "Preconditioner Block Symmetric Gauss Seidel does not support this matrix" );
            this->print ( );
            exit ( -1 );
        }

        template <typename ValueType>
        void lMatrix<ValueType>::BlocksPsgauss_seidel ( const lVector<ValueType> &invec, lVector<ValueType> *outvec,
                                                        const int num_blocks ) const
        {
            LOG_ERROR ( "Preconditioner Blocks Symmetric Gauss Seidel does not support this matrix" );
            this->print ( );
            exit ( -1 );
        }

        template <typename ValueType>
        void lMatrix<ValueType>::ZeroRows ( const int *index_set,
                                            const int size,
                                            const ValueType alpha )
        {
            LOG_ERROR ( "lMatrix::ZeroRows() does not support this matrix" );
            this->print ( );
            exit ( -1 );
        }

        template <typename ValueType>
        void lMatrix<ValueType>::ZeroCols ( const int *index_set,
                                            const int size,
                                            const ValueType alpha )
        {
            LOG_ERROR ( "lMatrix::ZeroCols() does not support this matrix" );
            this->print ( );
            exit ( -1 );
        }

        template <typename ValueType>
        void lMatrix<ValueType>::init_structure ( const int *rows, const int *cols )
        {
            LOG_ERROR ( "lMatrix::init_structure() does not support this matrix" );
            this->print ( );
            exit ( -1 );
        }

        template <typename ValueType>
        void lMatrix<ValueType>::add_value ( const int row, const int col, const ValueType val )
        {
            LOG_ERROR ( "lMatrix::add_value() does not support this matrix" );
            this->print ( );
            exit ( -1 );
        }

        template <typename ValueType>
        void lMatrix<ValueType>::add_values ( const int* rows, int num_rows, const int* cols, int num_cols, const ValueType* values )
        {
            LOG_ERROR ( "lMatrix::add_value() does not support this matrix" );
            this->print ( );
            exit ( -1 );
        }

        template <typename ValueType>
        void lMatrix<ValueType>::get_value ( const int row, const int col, ValueType *val ) const
        {
            LOG_ERROR ( "lMatrix::get_value() does not support this matrix" );
            this->print ( );
            exit ( -1 );
        }

        template <typename ValueType>
        void lMatrix<ValueType>::get_add_values ( const int* rows,
                                                  int num_rows,
                                                  const int* cols,
                                                  int num_cols,
                                                  const int* cols_target,
                                                  int num_cols_target,
                                                  ValueType* values ) const
        {
            LOG_ERROR ( "lMatrix::get_add_values() does not support this matrix" );
            this->print ( );
            exit ( -1 );
        }

        template <typename ValueType>
        void lMatrix<ValueType>::VectorMultAdd_submatrix ( const int* rows,
                                                           int num_rows,
                                                           const int* cols,
                                                           int num_cols,
                                                           const int* cols_input,
                                                           const ValueType* in_values,
                                                           ValueType* out_values ) const
        {
            LOG_ERROR ( "lMatrix::VectorMultAdd_submatrix() does not support this matrix" );
            this->print ( );
            exit ( -1 );
        }

        template <typename ValueType>
        void lMatrix<ValueType>::VectorMultAdd_submatrix_vanka ( const int* rows,
                                                                 int num_rows,
                                                                 const hiflow::la::lVector< ValueType > &invec,
                                                                 ValueType* out_values ) const
        {
            LOG_ERROR ( "lMatrix::VectorMultAdd_submatrix_vanka() does not support this matrix" );
            this->print ( );
            exit ( -1 );
        }

        template <typename ValueType>
        lMatrix<ValueType> *lMatrix<ValueType>::extract_submatrix ( const int start_row, const int start_col,
                                                                    const int end_row, const int end_col ) const
        {
            lMatrix<ValueType> *cpu_matrix = new CPUsimple_CSR_lMatrix<ValueType>;
            lMatrix<ValueType> *cpu_sub_matrix;
            lMatrix<ValueType> *sub_matrix;

            cpu_matrix->CloneFrom ( *this );

            cpu_sub_matrix = cpu_matrix->extract_submatrix ( start_row, start_col,
                                                             end_row, end_col );

            std::string sub_mat_name;
            sub_mat_name = "sub matrix from ";
            sub_mat_name.append ( this->get_name ( ) );

            sub_matrix = this->CloneWithoutContent ( );
            sub_matrix->Init ( cpu_sub_matrix->get_nnz ( ),
                               cpu_sub_matrix->get_num_row ( ),
                               cpu_sub_matrix->get_num_col ( ),
                               sub_mat_name );

            //  sub_matrix = init_matrix<ValueType>(cpu_sub_matrix->get_nnz(),
            //                                      cpu_sub_matrix->get_num_row(),
            //                                      cpu_sub_matrix->get_num_col(),
            //                                      sub_mat_name,
            //                                      this->get_platform(),
            //                                      this->get_implementation(),
            //                                      this->get_matrix_format());

            sub_matrix->CloneFrom ( *cpu_sub_matrix );

            delete cpu_matrix;
            delete cpu_sub_matrix;

            return sub_matrix;
        }

        template <typename ValueType>
        lMatrix<ValueType> * lMatrix<ValueType>::MatrixMult ( const lMatrix<ValueType> &inmat ) const
        {
            LOG_ERROR ( "lMatrix::MatrixMult() does not support this matrix" );
            this->print ( );
            exit ( -1 );
            return NULL;
        }

        template <typename ValueType>
        lMatrix<ValueType> * lMatrix<ValueType>::MatrixMultSupStructure ( const lMatrix<ValueType> &inmat ) const
        {
            LOG_ERROR ( "lMatrix::MatrixMultSupStructure() does not support this matrix" );
            this->print ( );
            exit ( -1 );
            return NULL;
        }

        template <typename ValueType>
        hiflow::la::lMatrix<ValueType> *lMatrix<ValueType>::MatrixSupSPower ( const int p ) const
        {
            LOG_ERROR ( "ERROR lMatrix::MatrixSupSPower() does not support this matrix" );
            this->print ( );
            exit ( -1 );
            return NULL;
        }

        template <typename ValueType>
        void lMatrix<ValueType>::MatrixAdd ( const lMatrix<ValueType> &inmat )
        {
            LOG_ERROR ( "lMatrix::MatrixAdd() does not support this matrix" );
            this->print ( );
            exit ( -1 );
        }

        template <typename ValueType>
        void lMatrix<ValueType>::Scale ( const ValueType alpha )
        {
            LOG_ERROR ( "lMatrix::Scale() does not support this matrix" );
            this->print ( );
            exit ( -1 );
        }

        template <typename ValueType>
        void lMatrix<ValueType>::ScaleOffdiag ( const ValueType alpha )
        {
            LOG_ERROR ( "lMatrix::ScaleOffdiag() does not support this matrix" );
            this->print ( );
            exit ( -1 );
        }

        template <typename ValueType>
        void lMatrix<ValueType>::GershgorinSpectrum ( ValueType *lambda_min, ValueType *lambda_max ) const
        {
            LOG_ERROR ( "lMatrix::GershroinSpectrum() does not support this matrix" );
            this->print ( );
            exit ( -1 );
        }

        template <typename ValueType>
        void lMatrix<ValueType>::Reorder ( const int *index )
        {
            LOG_ERROR ( "lMatrix::Reorder() does not support this matrix" );
            this->print ( );
            exit ( -1 );
        }

        template <typename ValueType>
        void lMatrix<ValueType>::extract_diagelements ( const int start_i, const int end_i, lVector<ValueType> *vec ) const
        {
            LOG_ERROR ( "lMatrix::extract_diagelements() does not support this matrix" );
            this->print ( );
            exit ( -1 );
        }

        template <typename ValueType>
        void lMatrix<ValueType>::extract_invdiagelements ( const int start_i, const int end_i, lVector<ValueType> *vec ) const
        {
            LOG_ERROR ( "lMatrix::extract_invdiagelements() does not support this matrix" );
            this->print ( );
            exit ( -1 );
        }

        template <typename ValueType>
        void lMatrix<ValueType>::Levelscheduling ( int &nlevels, int **level_sizes, int **permut_index ) const
        {
            LOG_ERROR ( "lMatrix::Levelscheduling() does not support this matrix" );
            this->print ( );
            exit ( -1 );
        }

        template <typename ValueType>
        void lMatrix<ValueType>::Multicoloring ( int &ncolors, int **color_sizes, int **permut_index ) const
        {
            LOG_ERROR ( "lMatrix::Multicoloring() does not support this matrix" );
            this->print ( );
            exit ( -1 );
        }

        template <typename ValueType>
        void lMatrix<ValueType>::Residual ( const lVector<ValueType> &b,
                                            const lVector<ValueType> &x,
                                            lVector<ValueType> *res ) const
        {

            this->VectorMult ( x, res );
            res->ScaleAdd ( -1.0, b );

        }

        template <typename ValueType>
        void lMatrix<ValueType>::CloneFrom ( const lMatrix<ValueType> &mat2 )
        {
            if ( this != &mat2 )
            {
                // if it is not empty init() will clean it
                this->Init ( mat2.get_nnz ( ), mat2.get_num_row ( ), mat2.get_num_col ( ), mat2.get_name ( ) );

                this->CopyStructureFrom ( mat2 );

                this->CopyFrom ( mat2 );
            }
        }

        template <typename ValueType>
        void lMatrix<ValueType>::ilu0 ( void )
        {
            LOG_ERROR ( "lMatrix::ilu0() does not support this matrix" );
            this->print ( );
            exit ( -1 );
        }

        template <typename ValueType>
        void lMatrix<ValueType>::ilup ( const int p )
        {
            LOG_ERROR ( "lMatrix::ilup(int) does not support this matrix" );
            this->print ( );
            exit ( -1 );
        }

        template <typename ValueType>
        void lMatrix<ValueType>::ilup ( const lMatrix<ValueType> &mat, const int p )
        {
            LOG_ERROR ( "lMatrix::ilup(lMatrix) does not support this matrix" );
            this->print ( );
            exit ( -1 );
        }

        template <typename ValueType>
        void lMatrix<ValueType>::ilu_solve ( const lVector<ValueType> &invec,
                                             lVector<ValueType> *outvec ) const
        {
            LOG_ERROR ( "lMatrix::ilu_solve() does not support this matrix" );
            this->print ( );
            exit ( -1 );
        }

        template <typename ValueType>
        void lMatrix<ValueType>::compress_me ( )
        {
            LOG_ERROR ( "ERROR lMatrix::compress_me() does not support this matrix" );
            this->print ( );
            exit ( -1 );
        }

        template <typename ValueType>
        bool lMatrix<ValueType>::issymmetric ( )
        {
            LOG_ERROR ( "ERROR lMatrix::issymmetric() does not support this matrix" );
            this->print ( );
            exit ( -1 );
            return false;
        }

        template <typename ValueType>
        void lMatrix<ValueType>::delete_diagonal ( )
        {
            LOG_ERROR ( "ERROR lMatrix::delete_diagonal() does not support this matrix" );
            this->print ( );
            exit ( -1 );
        }

        template <typename ValueType>
        void lMatrix<ValueType>::delete_offdiagonal ( )
        {
            LOG_ERROR ( "ERROR lMatrix::delete_offdiagonal() does not support this matrix" );
            this->print ( );
            exit ( -1 );
        }

        template <typename ValueType>
        void lMatrix<ValueType>::delete_lower_triangular ( )
        {
            LOG_ERROR ( "ERROR lMatrix::delete_lower_triangular() does not support this matrix" );
            this->print ( );
            exit ( -1 );
        }

        template <typename ValueType>
        void lMatrix<ValueType>::delete_strictly_lower_triangular ( )
        {
            LOG_ERROR ( "ERROR lMatrix::delete_strictly_lower_triangular() does not support this matrix" );
            this->print ( );
            exit ( -1 );
        }

        template <typename ValueType>
        void lMatrix<ValueType>::delete_upper_triangular ( )
        {
            LOG_ERROR ( "ERROR lMatrix::delete_upper_triangular() does not support this matrix" );
            this->print ( );
            exit ( -1 );
        }

        template <typename ValueType>
        void lMatrix<ValueType>::delete_strictly_upper_triangular ( )
        {
            LOG_ERROR ( "ERROR lMatrix::delete_strictly_upper_triangular() does not support this matrix" );
            this->print ( );
            exit ( -1 );
        }

        template <typename ValueType>
        void lMatrix<ValueType>::transpose_me ( )
        {
            LOG_ERROR ( "ERROR lMatrix::transpose_me() does not support this matrix" );
            this->print ( );
            exit ( -1 );
        }

        template class lMatrix<float>;
        template class lMatrix<double>;

    } // namespace la
} // namespace hiflow
