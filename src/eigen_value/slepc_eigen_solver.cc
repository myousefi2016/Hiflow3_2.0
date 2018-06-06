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

/// @author Philipp Gerstner

#include "slepc_eigen_solver.h"
#include "slepc.h"
#include "../linear_algebra/petsc_wrapper.h"
#include <vector>

namespace hiflow
{
    namespace la
    {

        void Complex2Real ( const std::vector<PetscScalar>& vec, std::vector<PetscReal>& real )
        {
            real.resize ( vec.size ( ), 0. );
            for ( int l = 0; l < vec.size ( ); ++l )
                real[l] = PetscRealPart ( vec[l] );
        }

        void Complex2Imag ( const std::vector<PetscScalar>& vec, std::vector<PetscReal>& imag )
        {
            imag.resize ( vec.size ( ), 0. );
            for ( int l = 0; l < vec.size ( ); ++l )
                imag[l] = PetscImaginaryPart ( vec[l] );
        }

        void Real2Complex ( const std::vector<PetscReal>& real, std::vector<PetscScalar>& vec, double s )
        {
            if ( s != 0 )
            {
                assert ( vec.size ( ) == real.size ( ) );
            }
            else
            {
                vec.resize ( real.size ( ), 0. );
                s = 1.;
            }
            for ( int l = 0; l < vec.size ( ); ++l )
            {
                vec[l] += s * real[l];
            }
        }

#ifdef WITH_COMPLEX_PETSC

        void Imag2Complex ( const std::vector<PetscReal>& imag, std::vector<PetscScalar>& vec, double s )
        {
            if ( s != 0 )
            {
                assert ( vec.size ( ) == imag.size ( ) );
            }
            else
            {
                vec.resize ( imag.size ( ), 0. );
                s = 1.;
            }
            for ( int l = 0; l < vec.size ( ); ++l )
            {
                vec[l] += s * PETSC_i * imag[l];
            }
        }
#endif

        /// Matrix Vector product of (Hiflow) matrix with PETSc vector

        template <class LAD>
        PetscErrorCode MatVecMult ( Mat A, Vec x, Vec y )
        {
            void *ctx;
            PetscErrorCode ierr = 0;
            Operator_wrapper<LAD> op_ptr;
            std::vector<int> ind;
            std::vector<typename LAD::DataType> vals;
            std::vector<PetscScalar> petsc_x;
            std::vector<PetscScalar> petsc_y;
            int size;

            // Get required data structures	
            ierr = MatShellGetContext ( A, &ctx );
            if ( ierr != 0 ) return ierr;

            op_ptr = *( Operator_wrapper<LAD>* ) ctx;

            int dof_offset = op_ptr.dof_offset_;
            int local_dofs = op_ptr.local_dofs_;

            ierr = VecGetLocalSize ( x, &size );
            if ( ierr != 0 ) return ierr;
            assert ( size == local_dofs );

            // PetscVec x_r -> HiflowVec vec_ 
            petsc_x.resize ( size, 0. );
            ind.resize ( size, 0 );
            op_ptr.vec_->Zeros ( );
            for ( int l = 0; l < size; l++ )
            {
                ind[l] = dof_offset + l;
            }

            ierr = VecGetValues ( x, size, &ind[0], &petsc_x[0] );
            if ( ierr != 0 ) return ierr;

#ifdef WITH_COMPLEX_PETSC 
            Complex2Real ( petsc_x, vals );
            op_ptr.vec_->SetValues ( &ind[0], size, &vals[0] );
#else 
            op_ptr.vec_->SetValues ( &ind[0], size, &petsc_x[0] );
#endif

            op_ptr.vec_->UpdateGhost ( );

            // HiflowVec res = HiflowOp * HiflowVec vec_ 	
            typename LAD::VectorType res;
            res.CloneFromWithoutContent ( *op_ptr.vec_ );
            op_ptr.op_->VectorMult ( *op_ptr.vec_, &res );

            // HiflowVec res -> PetscVec y
            // y = A_r * x_r	
            vals.resize ( size, 0. );
            res.GetValues ( &ind[0], size, &vals[0] );

#ifdef WITH_COMPLEX_PETSC 
            Real2Complex ( vals, petsc_y, 0. );
#else  
            ierr = VecSetValues ( y, size, &ind[0], &vals[0], INSERT_VALUES );
            if ( ierr != 0 ) return ierr;
            ierr = VecAssemblyBegin ( y );
            if ( ierr != 0 ) return ierr;
            ierr = VecAssemblyEnd ( y );
            if ( ierr != 0 ) return ierr;
#endif

#ifdef WITH_COMPLEX_PETSC
            if ( op_ptr.opi_ != NULL )
            {
                // HiflowVec res = HiflowOp_i * HiflowVec vec_
                res.CloneFromWithoutContent ( *op_ptr.vec_ );
                op_ptr.opi_->VectorMult ( *op_ptr.vec_, &res );

                // y = A_r * x_r + i * A_i * x_r	
                vals.resize ( size, 0. );
                res.GetValues ( &ind[0], size, &vals[0] );
                Imag2Complex ( vals, petsc_y, 1 );
            }

            // PetscVec x_i -> HiflowVec vec_ 
            Complex2Imag ( petsc_x, vals );

            op_ptr.vec_->Zeros ( );
            op_ptr.vec_->SetValues ( &ind[0], size, &vals[0] );
            op_ptr.vec_->UpdateGhost ( );

            // HiflowVec res = HiflowOp * HiflowVec vec_ 	
            res.CloneFromWithoutContent ( *op_ptr.vec_ );
            op_ptr.op_->VectorMult ( *op_ptr.vec_, &res );

            // HiflowVec res -> PetscVec y
            // y = A_r * x_r + i * (A_i * x_r + A_r * x_i)	
            vals.resize ( size, 0. );
            res.GetValues ( &ind[0], size, &vals[0] );
            Imag2Complex ( vals, petsc_y, 1. );

            if ( op_ptr.opi_ != NULL )
            {
                // HiflowVec res = HiflowOp_i * HiflowVec vec_
                res.CloneFromWithoutContent ( *op_ptr.vec_ );
                op_ptr.opi_->VectorMult ( *op_ptr.vec_, &res );

                // y = A_r * x_r + - A_i * x_i +  i * (A_i * x_i + A_r * x_i)	
                vals.resize ( size, 0. );
                res.GetValues ( &ind[0], size, &vals[0] );
                Real2Complex ( vals, petsc_y, -1. );

                ierr = VecSetValues ( y, size, &ind[0], &petsc_y[0], INSERT_VALUES );
                if ( ierr != 0 ) return ierr;
            }

            ierr = VecAssemblyBegin ( y );
            if ( ierr != 0 ) return ierr;
            ierr = VecAssemblyEnd ( y );
            if ( ierr != 0 ) return ierr;
#endif  
            return ierr;
        }

        /// Matrix Vector product of transposed (Hiflow) matrix with PETSc vector

        template <class LAD>
        PetscErrorCode MatTransVecMult ( Mat A, Vec x, Vec y )
        {
            void *ctx;
            PetscErrorCode ierr = 0;
            Operator_wrapper<LAD> op_ptr;
            std::vector<int> ind;
            std::vector<typename LAD::DataType> vals;
            std::vector<PetscScalar> petsc_x;
            std::vector<PetscScalar> petsc_y;
            int size;

            // Get required data structures	
            ierr = MatShellGetContext ( A, &ctx );
            if ( ierr != 0 ) return ierr;

            op_ptr = *( Operator_wrapper<LAD>* ) ctx;

            int dof_offset = op_ptr.dof_offset_;
            int local_dofs = op_ptr.local_dofs_;

            ierr = VecGetLocalSize ( x, &size );
            if ( ierr != 0 ) return ierr;
            assert ( size == local_dofs );

            // PetscVec x_r -> HiflowVec vec_ 
            petsc_x.resize ( size, 0. );
            ind.resize ( size, 0 );
            op_ptr.vec_->Zeros ( );
            for ( int l = 0; l < size; l++ )
            {
                ind[l] = dof_offset + l;
            }

            ierr = VecGetValues ( x, size, &ind[0], &petsc_x[0] );
            if ( ierr != 0 ) return ierr;

#ifdef WITH_COMPLEX_PETSC 
            Complex2Real ( petsc_x, vals );
            op_ptr.vec_->SetValues ( &ind[0], size, &vals[0] );
#else 
            op_ptr.vec_->SetValues ( &ind[0], size, &petsc_x[0] );
#endif

            op_ptr.vec_->UpdateGhost ( );

            // HiflowVec res = HiflowOp * HiflowVec vec_ 	
            typename LAD::VectorType res;
            res.CloneFromWithoutContent ( *op_ptr.vec_ );
            op_ptr.opT_->VectorMult ( *op_ptr.vec_, &res );

            // HiflowVec res -> PetscVec y
            // y = A_r * x_r	
            vals.resize ( size, 0. );
            res.GetValues ( &ind[0], size, &vals[0] );

#ifdef WITH_COMPLEX_PETSC 
            Real2Complex ( vals, petsc_y, 0. );
#else  
            ierr = VecSetValues ( y, size, &ind[0], &vals[0], INSERT_VALUES );
            if ( ierr != 0 ) return ierr;
            ierr = VecAssemblyBegin ( y );
            if ( ierr != 0 ) return ierr;
            ierr = VecAssemblyEnd ( y );
            if ( ierr != 0 ) return ierr;
#endif

#ifdef WITH_COMPLEX_PETSC
            if ( op_ptr.opiT_ != NULL )
            {
                // HiflowVec res = HiflowOp_i * HiflowVec vec_
                res.CloneFromWithoutContent ( *op_ptr.vec_ );
                op_ptr.opiT_->VectorMult ( *op_ptr.vec_, &res );

                // y = A_r * x_r + i * A_i * x_r	
                vals.resize ( size, 0. );
                res.GetValues ( &ind[0], size, &vals[0] );
                Imag2Complex ( vals, petsc_y, 1 );
            }

            // PetscVec x_i -> HiflowVec vec_ 
            Complex2Imag ( petsc_x, vals );

            op_ptr.vec_->Zeros ( );
            op_ptr.vec_->SetValues ( &ind[0], size, &vals[0] );
            op_ptr.vec_->UpdateGhost ( );

            // HiflowVec res = HiflowOp * HiflowVec vec_ 	
            res.CloneFromWithoutContent ( *op_ptr.vec_ );
            op_ptr.opT_->VectorMult ( *op_ptr.vec_, &res );

            // HiflowVec res -> PetscVec y
            // y = A_r * x_r + i * (A_i * x_r + A_r * x_i)	
            vals.resize ( size, 0. );
            res.GetValues ( &ind[0], size, &vals[0] );
            Imag2Complex ( vals, petsc_y, 1. );

            if ( op_ptr.opiT_ != NULL )
            {
                // HiflowVec res = HiflowOp_i * HiflowVec vec_
                res.CloneFromWithoutContent ( *op_ptr.vec_ );
                op_ptr.opiT_->VectorMult ( *op_ptr.vec_, &res );

                // y = A_r * x_r + - A_i * x_i +  i * (A_i * x_i + A_r * x_i)	
                vals.resize ( size, 0. );
                res.GetValues ( &ind[0], size, &vals[0] );
                Real2Complex ( vals, petsc_y, -1. );

                ierr = VecSetValues ( y, size, &ind[0], &petsc_y[0], INSERT_VALUES );
                if ( ierr != 0 ) return ierr;
            }

            ierr = VecAssemblyBegin ( y );
            if ( ierr != 0 ) return ierr;
            ierr = VecAssemblyEnd ( y );
            if ( ierr != 0 ) return ierr;
#endif  
            return ierr;
        }

        namespace slepc
        {

            struct EPS_wrapper
            {
                ::EPS eps_;
            };

            struct ST_wrapper
            {
                ::ST st_;
            };

            struct SVD_wrapper
            {
                ::SVD svd_;
            };

            struct KSP_wrapper
            {
                ::KSP ksp_;
            };
        };

        template <class LAD>
        SLEPcEigenSolver<LAD>::SLEPcEigenSolver ( slepc::EPS_wrapper* ptr_eps_wrapper )
        : ptr_eps_wrapper_ ( ptr_eps_wrapper ),
        comm_ ( MPI_COMM_NULL ),
        rank_ ( -1 ),
        nb_procs_ ( -1 ),
        mat_A_set_ ( false ),
        mat_B_set_ ( false ),
        ptr_A_wrapper_ ( new petsc::Mat_wrapper ),
        ptr_B_wrapper_ ( new petsc::Mat_wrapper ),
        EigenSolver<LAD>( )
        {
            struct_A_.is_set = false;
            struct_Ai_.is_set = false;
            struct_B_.is_set = false;
            struct_Bi_.is_set = false;
        }

        template <class LAD>
        SLEPcEigenSolver<LAD>::SLEPcEigenSolver ( slepc::SVD_wrapper* ptr_svd_wrapper )
        : ptr_svd_wrapper_ ( ptr_svd_wrapper ),
        rank_ ( -1 ),
        nb_procs_ ( -1 ),
        comm_ ( MPI_COMM_NULL ),
        mat_A_set_ ( false ),
        mat_B_set_ ( false ),
        ptr_A_wrapper_ ( new petsc::Mat_wrapper ),
        ptr_B_wrapper_ ( new petsc::Mat_wrapper ),
        EigenSolver<LAD>( )
        {
            struct_A_.is_set = false;
            struct_Ai_.is_set = false;
            struct_B_.is_set = false;
            struct_Bi_.is_set = false;
        }

        template <class LAD>
        void SLEPcEigenSolver<LAD>::SetupMatrixStructure ( const Couplings<DataType>& cp, const std::vector<int>& d_rows, const std::vector<int>& d_cols,
                                                           const std::vector<int>& o_rows, const std::vector<int>& o_cols, slepc::MatrixStructure& m_struct )
        {
            assert ( m_struct.mat_type == petsc::MPI_SPARSE );

            m_struct.diag_cols.resize ( d_cols.size ( ), 0. );
            m_struct.diag_rows.resize ( d_rows.size ( ), 0. );
            m_struct.off_diag_cols.resize ( o_cols.size ( ), 0. );
            m_struct.off_diag_rows.resize ( o_rows.size ( ), 0. );

            for ( int l = 0; l < d_cols.size ( ); ++l )
                m_struct.diag_cols[l] = d_cols[l];

            for ( int l = 0; l < d_rows.size ( ); ++l )
                m_struct.diag_rows[l] = d_rows[l];

            for ( int l = 0; l < o_cols.size ( ); ++l )
                m_struct.off_diag_cols[l] = o_cols[l];

            for ( int l = 0; l < o_rows.size ( ); ++l )
                m_struct.off_diag_rows[l] = o_rows[l];

            m_struct.diag_nnz = m_struct.diag_cols.size ( );
            m_struct.off_diag_nnz = m_struct.off_diag_cols.size ( );

            m_struct.local_dofs = cp.nb_dofs ( this->rank_ );
            m_struct.dof_offset = cp.dof_offset ( this->rank_ );
            m_struct.d_nnz.resize ( m_struct.local_dofs, 0 );
            m_struct.o_nnz.resize ( m_struct.local_dofs, 0 );

            for ( int l = 0; l < m_struct.diag_nnz; ++l )
                m_struct.d_nnz[m_struct.diag_rows[l] - m_struct.dof_offset] += 1;

            for ( int l = 0; l < m_struct.off_diag_nnz; ++l )
                m_struct.o_nnz[m_struct.off_diag_rows[l] - m_struct.dof_offset] += 1;

            m_struct.nz_col_ind_in_diag_row.resize ( m_struct.local_dofs );
            m_struct.nz_col_ind_in_offdiag_row.resize ( m_struct.local_dofs );

            for ( int l = 0; l < m_struct.diag_nnz; ++l )
                m_struct.nz_col_ind_in_diag_row[m_struct.diag_rows[l] - m_struct.dof_offset].push_back ( m_struct.diag_cols[l] );

            for ( int l = 0; l < m_struct.off_diag_nnz; ++l )
                m_struct.nz_col_ind_in_offdiag_row[m_struct.off_diag_rows[l] - m_struct.dof_offset].push_back ( m_struct.off_diag_cols[l] );

            m_struct.is_set = true;
        }

        template <class LAD>
        void SLEPcEigenSolver<LAD>::SetupMatrixStructure ( const Couplings<DataType>& cp, slepc::MatrixStructure& m_struct )
        {
            assert ( m_struct.mat_type == petsc::MPI_MATFREE );

            m_struct.local_dofs = cp.nb_dofs ( this->rank_ );
            m_struct.dof_offset = cp.dof_offset ( this->rank_ );
            m_struct.is_set = true;
            InitVec ( cp );
        }

        template <class LAD>
        void SLEPcEigenSolver<LAD>::ClearMatrixStructure ( slepc::MatrixStructure& m_struct )
        {
            m_struct.diag_cols.clear ( );
            m_struct.diag_rows.clear ( );
            m_struct.off_diag_cols.clear ( );
            m_struct.off_diag_rows.clear ( );
            m_struct.d_nnz.clear ( );
            m_struct.o_nnz.clear ( );

            for ( int j = 0; j < m_struct.nz_col_ind_in_diag_row.size ( ); ++j )
            {
                m_struct.nz_col_ind_in_diag_row[j].clear ( );
                m_struct.nz_col_ind_in_offdiag_row[j].clear ( );
            }
            m_struct.nz_col_ind_in_diag_row.clear ( );
            m_struct.nz_col_ind_in_offdiag_row.clear ( );

            m_struct.is_set = false;
        }

        template <>
        void SLEPcEigenSolver<LADescriptorCoupledD>::InitVec ( const Couplings<DataType>& cp )
        {
            aux_vec_.Init ( this->comm_, cp, CPU, NAIVE );
            aux_vec_.InitStructure ( );
        }

        template <>
        void SLEPcEigenSolver<LADescriptorHypreD>::InitVec ( const Couplings<DataType>& cp )
        {
            aux_vec_.Init ( this->comm_, cp );
        }

#ifndef WITH_COMPLEX_PETSC

        template <>
        void SLEPcEigenSolver<LADescriptorPETSc>::InitVec ( const Couplings<DataType>& cp )
        {
            aux_vec_.Init ( this->comm_, cp );
        }
#endif

        template <class LAD>
        void SLEPcEigenSolver<LAD>::SetupOperatorA ( const OperatorType* opA, const OperatorType* opAi,
                                                     const Couplings<DataType>& cp, const std::vector<int>& d_rows, const std::vector<int>& d_cols, const std::vector<int>& o_rows, const std::vector<int>& o_cols )
        {
            struct_A_.mat_type = petsc::MPI_SPARSE;
            SetupMatrixStructure ( cp, d_rows, d_cols, o_rows, o_cols, struct_A_ );

            if ( opAi != NULL )
            {
                struct_Ai_.mat_type = petsc::MPI_SPARSE;
                SetupMatrixStructure ( cp, d_rows, d_cols, o_rows, o_cols, struct_Ai_ );
                this->complex_A_ = true;
            }
            SetupPetscMat ( opA, opAi, struct_A_, struct_Ai_, ptr_A_wrapper_ );
            mat_A_set_ = true;
        }

        template <class LAD>
        void SLEPcEigenSolver<LAD>::SetupOperatorA ( const OperatorType* opA, const OperatorType* opAt, const OperatorType* opAi, const OperatorType* opAit, const Couplings<DataType>& cp )
        {
            struct_A_.mat_type = petsc::MPI_MATFREE;
            SetupMatrixStructure ( cp, struct_A_ );

            if ( opAi != NULL )
            {
                struct_Ai_.mat_type = petsc::MPI_MATFREE;
                SetupMatrixStructure ( cp, struct_Ai_ );
                this->complex_A_ = true;
            }
            SetupPetscMatShell ( opA, opAt, opAi, opAit, ptr_A_wrapper_, 0 );
            mat_A_set_ = true;
        }

        template <class LAD>
        void SLEPcEigenSolver<LAD>::SetupOperatorB ( const OperatorType* opB, const OperatorType* opBi,
                                                     const Couplings<DataType>& cp, const std::vector<int>& d_rows, const std::vector<int>& d_cols, const std::vector<int>& o_rows, const std::vector<int>& o_cols )
        {
            struct_B_.mat_type = petsc::MPI_SPARSE;
            SetupMatrixStructure ( cp, d_rows, d_cols, o_rows, o_cols, struct_B_ );

            if ( opBi != NULL )
            {
                struct_Bi_.mat_type = petsc::MPI_SPARSE;
                SetupMatrixStructure ( cp, d_rows, d_cols, o_rows, o_cols, struct_Bi_ );
                this->complex_B_ = true;
            }

            //	std::cout << opB << " " << opBi << " " << struct_B_.is_set << " " << struct_Bi_.is_set << std::endl;

            SetupPetscMat ( opB, opBi, struct_B_, struct_Bi_, ptr_B_wrapper_ );
            mat_B_set_ = true;
        }

        template <class LAD>
        void SLEPcEigenSolver<LAD>::SetupOperatorB ( const OperatorType* opB, const OperatorType* opBt, const OperatorType* opBi, const OperatorType* opBit, const Couplings<DataType>& cp )
        {
            struct_B_.mat_type = petsc::MPI_MATFREE;
            SetupMatrixStructure ( cp, struct_B_ );

            if ( opBi != NULL )
            {
                struct_Bi_.mat_type = petsc::MPI_MATFREE;
                SetupMatrixStructure ( cp, struct_Bi_ );
                this->complex_B_ = true;
            }
            SetupPetscMatShell ( opB, opBt, opBi, opBit, ptr_B_wrapper_, 1 );
            mat_B_set_ = true;
        }

        template <class LAD>
        void SLEPcEigenSolver<LAD>::SetupPetscMat ( const OperatorType* op, const OperatorType* opi, slepc::MatrixStructure struct_real, slepc::MatrixStructure struct_imag, petsc::Mat_wrapper* ptr_mat_wrapper )
        {

            // Create parallel sparse matrix	
            std::vector<int> d_nnz;
            std::vector<int> o_nnz;
            d_nnz.resize ( struct_real.d_nnz.size ( ), 0 );
            o_nnz.resize ( struct_real.o_nnz.size ( ), 0 );
            for ( int i = 0; i < d_nnz.size ( ); ++i )
            {
                if ( struct_real.is_set )
                    d_nnz[i] += struct_real.d_nnz[i];
                if ( struct_imag.is_set )
                    d_nnz[i] += struct_imag.d_nnz[i];
            }
            for ( int i = 0; i < o_nnz.size ( ); ++i )
            {
                if ( struct_real.is_set )
                    o_nnz[i] += struct_real.o_nnz[i];
                if ( struct_imag.is_set )
                    o_nnz[i] += struct_imag.o_nnz[i];
            }

            //std::cout << this->rank_ << ": " << local_dofs_ << 	" - " << op.num_rows_local() << " / " << op.num_cols_local() << std::endl;
            MatCreate ( this->comm_, &ptr_mat_wrapper->mat_ );
            MatSetType ( ptr_mat_wrapper->mat_, MATMPIAIJ );
            MatSetSizes ( ptr_mat_wrapper->mat_, op->num_rows_local ( ), op->num_cols_local ( ), op->num_rows_global ( ), op->num_cols_global ( ) );
            MatMPIAIJSetPreallocation ( ptr_mat_wrapper->mat_, 0, vec2ptr ( d_nnz ), 0, vec2ptr ( o_nnz ) );

            // Entry by entry copy	
            /*
            for (int l=0; l<diag_nnz_; ++l)
            {
                    DataType val;
                    op->GetValues(&diag_rows_[l], 1, &diag_cols_[l], 1, &val);
                    // std::cout << diag_rows_[l] << " , " << diag_cols_[l] << " : " << val << std::endl;			
                    MatSetValue(ptr_mat_wrapper->mat_, diag_rows_[l], diag_cols_[l], val, INSERT_VALUES);  
            }
            for (int l=0; l<off_diag_nnz_; ++l)
            {
                    DataType val;
                    op->GetValues(&off_diag_rows_[l], 1, &off_diag_cols_[l], 1, &val);
                    //std::cout << off_diag_rows_[l] << " , " << off_diag_cols_[l] << " : " << val << std::endl;				
                    MatSetValue(ptr_mat_wrapper->mat_, off_diag_rows_[l], off_diag_cols_[l], val, INSERT_VALUES);
            }
             */
            // Row by row copy of real part
            if ( struct_real.is_set )
            {
                for ( int j = 0; j < struct_real.local_dofs; ++j )
                {
                    std::vector<DataType> diag_vals;
                    std::vector<DataType> offdiag_vals;

                    std::vector<PetscScalar> petsc_diag_vals;
                    std::vector<PetscScalar> petsc_offdiag_vals;

                    int nnz_diag_row = struct_real.d_nnz[j];
                    int nnz_offdiag_row = struct_real.o_nnz[j];
                    diag_vals.resize ( nnz_diag_row, 0. );
                    offdiag_vals.resize ( nnz_offdiag_row, 0. );
                    petsc_diag_vals.resize ( nnz_diag_row, 0. );
                    petsc_offdiag_vals.resize ( nnz_offdiag_row, 0. );

                    int cur_row = struct_real.dof_offset + j;

                    op->GetValues ( &cur_row, 1, &struct_real.nz_col_ind_in_diag_row[j][0], nnz_diag_row, &diag_vals[0] );
                    op->GetValues ( &cur_row, 1, &struct_real.nz_col_ind_in_offdiag_row[j][0], nnz_offdiag_row, &offdiag_vals[0] );

                    for ( int l = 0; l < nnz_diag_row; ++l )
                        petsc_diag_vals[l] = diag_vals[l];

                    for ( int l = 0; l < nnz_offdiag_row; ++l )
                        petsc_offdiag_vals[l] = offdiag_vals[l];

                    MatSetValues ( ptr_mat_wrapper->mat_, 1, &cur_row, nnz_diag_row, &struct_real.nz_col_ind_in_diag_row[j][0], &petsc_diag_vals[0], INSERT_VALUES );
                    MatSetValues ( ptr_mat_wrapper->mat_, 1, &cur_row, nnz_offdiag_row, &struct_real.nz_col_ind_in_offdiag_row[j][0], &petsc_offdiag_vals[0], INSERT_VALUES );
                }
                if ( struct_imag.is_set )
                {
                    MatAssemblyBegin ( ptr_mat_wrapper->mat_, MAT_FLUSH_ASSEMBLY );
                    MatAssemblyEnd ( ptr_mat_wrapper->mat_, MAT_FLUSH_ASSEMBLY );
                }
                else
                {
                    MatAssemblyBegin ( ptr_mat_wrapper->mat_, MAT_FINAL_ASSEMBLY );
                    MatAssemblyEnd ( ptr_mat_wrapper->mat_, MAT_FINAL_ASSEMBLY );
                }
            }

#ifdef WITH_COMPLEX_PETSC
            if ( struct_imag.is_set )
            {
                // Row by row copy of imag part part 
                for ( int j = 0; j < struct_imag.local_dofs; ++j )
                {
                    std::vector<DataType> diag_vals;
                    std::vector<DataType> offdiag_vals;
                    std::vector<PetscScalar> petsc_diag_vals;
                    std::vector<PetscScalar> petsc_offdiag_vals;
                    int nnz_diag_row = struct_imag.d_nnz[j];
                    int nnz_offdiag_row = struct_imag.o_nnz[j];
                    diag_vals.resize ( nnz_diag_row, 0. );
                    offdiag_vals.resize ( nnz_offdiag_row, 0. );

                    int cur_row = struct_imag.dof_offset + j;

                    opi->GetValues ( &cur_row, 1, &struct_imag.nz_col_ind_in_diag_row[j][0], nnz_diag_row, &diag_vals[0] );
                    opi->GetValues ( &cur_row, 1, &struct_imag.nz_col_ind_in_offdiag_row[j][0], nnz_offdiag_row, &offdiag_vals[0] );

                    petsc_diag_vals.resize ( diag_vals.size ( ), 0. );
                    petsc_offdiag_vals.resize ( offdiag_vals.size ( ), 0. );
                    for ( int l = 0; l < diag_vals.size ( ); ++l )
                    {
                        petsc_diag_vals[l] = PETSC_i * diag_vals[l];
                    }
                    for ( int l = 0; l < offdiag_vals.size ( ); ++l )
                    {
                        petsc_offdiag_vals[l] = PETSC_i * offdiag_vals[l];
                    }

                    MatSetValues ( ptr_mat_wrapper->mat_, 1, &cur_row, nnz_diag_row, &struct_imag.nz_col_ind_in_diag_row[j][0], &petsc_diag_vals[0], ADD_VALUES );
                    MatSetValues ( ptr_mat_wrapper->mat_, 1, &cur_row, nnz_offdiag_row, &struct_imag.nz_col_ind_in_offdiag_row[j][0], &petsc_offdiag_vals[0], ADD_VALUES );
                }
                MatAssemblyBegin ( ptr_mat_wrapper->mat_, MAT_FINAL_ASSEMBLY );
                MatAssemblyEnd ( ptr_mat_wrapper->mat_, MAT_FINAL_ASSEMBLY );
            }

#endif 
        }

        template <class LAD>
        void SLEPcEigenSolver<LAD>::SetupPetscMatShell ( const OperatorType* op, const OperatorType* opT, const OperatorType* opi, const OperatorType* opiT, petsc::Mat_wrapper* ptr_mat_wrapper, int matrix )
        {
            // Create operator that performs a matrix vector multiplication 
            if ( matrix == 0 )
            {
                opA_ptr_.op_ = op;
                opA_ptr_.opT_ = opT;
                opA_ptr_.opi_ = opi;
                opA_ptr_.opiT_ = opiT;
                opA_ptr_.local_dofs_ = struct_A_.local_dofs;
                opA_ptr_.dof_offset_ = struct_A_.dof_offset;
                opA_ptr_.vec_ = &aux_vec_;
                MatCreateShell ( this->comm_, op->num_rows_local ( ), op->num_cols_local ( ), op->num_rows_global ( ), op->num_cols_global ( ), &opA_ptr_, &ptr_mat_wrapper->mat_ );
            }
            else
            {
                opB_ptr_.op_ = op;
                opB_ptr_.opT_ = opT;
                opB_ptr_.opi_ = opi;
                opB_ptr_.opiT_ = opiT;
                opB_ptr_.local_dofs_ = struct_B_.local_dofs;
                opB_ptr_.dof_offset_ = struct_B_.dof_offset;
                opB_ptr_.vec_ = &aux_vec_;
                MatCreateShell ( this->comm_, op->num_rows_local ( ), op->num_cols_local ( ), op->num_rows_global ( ), op->num_cols_global ( ), &opB_ptr_, &ptr_mat_wrapper->mat_ );
            }
            MatShellSetOperation ( ptr_mat_wrapper->mat_, MATOP_MULT, ( void(* )( ) )MatVecMult<LAD> );
            MatShellSetOperation ( ptr_mat_wrapper->mat_, MATOP_MULT_TRANSPOSE, ( void(* )( ) )MatTransVecMult<LAD> );
        }

        template <class LAD>
        void SLEPcEigenSolver<LAD>::GetEigenValue ( int j, DataType& real, DataType& imag )
        {
            real = this->eig_val_[j].real_;
            imag = this->eig_val_[j].imag_;
        }

        template <class LAD>
        void SLEPcEigenSolver<LAD>::GetEigenVector ( int j, VectorType& v_real, VectorType& v_imag )
        {
            std::vector<int> ind;
            std::vector<PetscReal> r_vals;
            std::vector<PetscReal> i_vals;
            std::vector<PetscScalar> vals;
            int size;

#ifdef WITH_COMPLEX_PETSC
            VecAssemblyBegin ( this->ptr_v_wrapper_[j]->vec_ );
            VecAssemblyEnd ( this->ptr_v_wrapper_[j]->vec_ );

            VecGetLocalSize ( this->ptr_v_wrapper_[j]->vec_, &size );
            assert ( size == struct_A_.local_dofs );

            vals.resize ( size, 0. );
            r_vals.resize ( size, 0. );
            i_vals.resize ( size, 0. );
            ind.resize ( size, 0 );
            for ( int l = 0; l < size; l++ )
                ind[l] = struct_A_.dof_offset + l;

            VecGetValues ( this->ptr_v_wrapper_[j]->vec_, size, &ind[0], &vals[0] );
            for ( int j = 0; j < vals.size ( ); ++j )
            {
                r_vals[j] = PetscRealPart ( vals[j] );
                i_vals[j] = PetscImaginaryPart ( vals[j] );
            }

#else 
            VecAssemblyBegin ( this->ptr_vr_wrapper_[j]->vec_ );
            VecAssemblyBegin ( this->ptr_vi_wrapper_[j]->vec_ );
            VecAssemblyEnd ( this->ptr_vr_wrapper_[j]->vec_ );
            VecAssemblyEnd ( this->ptr_vi_wrapper_[j]->vec_ );

            VecGetLocalSize ( this->ptr_vr_wrapper_[j]->vec_, &size );
            assert ( size == struct_A_.local_dofs );

            r_vals.resize ( size, 0. );
            i_vals.resize ( size, 0. );
            ind.resize ( size, 0 );

            for ( int l = 0; l < size; l++ )
                ind[l] = struct_A_.dof_offset + l;

            VecGetValues ( this->ptr_vr_wrapper_[j]->vec_, size, &ind[0], &r_vals[0] );
            VecGetValues ( this->ptr_vi_wrapper_[j]->vec_, size, &ind[0], &i_vals[0] );
#endif

            v_real.SetValues ( &ind[0], size, &r_vals[0] );
            v_imag.SetValues ( &ind[0], size, &i_vals[0] );

            v_real.UpdateGhost ( );
            v_imag.UpdateGhost ( );

        }

        /*
        template <class LAD>
        void SLEPcEigenSolver<LAD>::GetSingularValue(int j, DataType& sig)
        {
                sig = this->sig_[j];
        }
         */

        /*
        template <class LAD>
        void SLEPcEigenSolver<LAD>::GetSingularVector(int j, VectorType& v_right, VectorType& v_left)
        {
                VecAssemblyBegin(this->ptr_vr_wrapper_[j]->vec_);
                VecAssemblyBegin(this->ptr_vi_wrapper_[j]->vec_);
                VecAssemblyEnd	(this->ptr_vr_wrapper_[j]->vec_);
                VecAssemblyEnd	(this->ptr_vi_wrapper_[j]->vec_);

                std::vector<int> ind;
                std::vector<DataType> r_vals;
                std::vector<DataType> i_vals;

                int size;
                VecGetLocalSize(this->ptr_vr_wrapper_[j]->vec_, &size);
                assert ( size == struct_A_.local_dofs );

                r_vals.resize(size, 0.);
                i_vals.resize(size, 0.);
                ind.resize(size, 0);

                for (int l=0; l<size;l++) 
                        ind[l] = struct_A_.dof_offset+l;

                VecGetValues(this->ptr_vr_wrapper_[j]->vec_, size, &ind[0], &r_vals[0]);
                VecGetValues(this->ptr_vi_wrapper_[j]->vec_, size, &ind[0], &i_vals[0]);	

                v_right.SetValues(&ind[0], size, &r_vals[0]);
                v_left.SetValues(&ind[0], size, &i_vals[0]);
        }
         */

        template <class LAD>
        void SLEPcEigenSolver<LAD>::Clear ( )
        {
            PetscErrorCode ierr;

            if ( mat_A_set_ )
                ierr = MatDestroy ( &ptr_A_wrapper_->mat_ );

            if ( mat_B_set_ )
                ierr = MatDestroy ( &ptr_B_wrapper_->mat_ );

            mat_A_set_ = false;
            mat_B_set_ = false;

            for ( int l = 0; l < ptr_v_wrapper_.size ( ); l++ )
                ierr = VecDestroy ( &ptr_v_wrapper_[l]->vec_ );

            for ( int l = 0; l < ptr_vr_wrapper_.size ( ); l++ )
                ierr = VecDestroy ( &ptr_vr_wrapper_[l]->vec_ );

            for ( int l = 0; l < ptr_vi_wrapper_.size ( ); l++ )
                ierr = VecDestroy ( &ptr_vi_wrapper_[l]->vec_ );

            eig_val_.clear ( );

            ptr_v_wrapper_.clear ( );
            ptr_vr_wrapper_.clear ( );
            ptr_vi_wrapper_.clear ( );

            ClearMatrixStructure ( struct_A_ );
            ClearMatrixStructure ( struct_Ai_ );
            ClearMatrixStructure ( struct_B_ );
            ClearMatrixStructure ( struct_Bi_ );

            if ( struct_A_.mat_type == petsc::MPI_MATFREE || struct_B_.mat_type == petsc::MPI_MATFREE )
                aux_vec_.Zeros ( );

            EigenSolver<LAD>::Clear ( );
        }

#ifdef WITH_PETSC  
        template class SLEPcEigenSolver<LADescriptorPETSc>;
        template PetscErrorCode MatVecMult<LADescriptorPETSc> ( Mat A, Vec x, Vec y );
        template PetscErrorCode MatTransVecMult<LADescriptorPETSc> ( Mat A, Vec x, Vec y );
#endif
#ifdef WITH_HYPRE
        template class SLEPcEigenSolver<LADescriptorHypreD>;
        template PetscErrorCode MatVecMult<LADescriptorHypreD> ( Mat A, Vec x, Vec y );
        template PetscErrorCode MatTransVecMult<LADescriptorHypreD> ( Mat A, Vec x, Vec y );
#endif
        template class SLEPcEigenSolver<LADescriptorCoupledD>;
        template PetscErrorCode MatVecMult<LADescriptorCoupledD> ( Mat A, Vec x, Vec y );
        template PetscErrorCode MatTransVecMult<LADescriptorCoupledD> ( Mat A, Vec x, Vec y );

    }
}
