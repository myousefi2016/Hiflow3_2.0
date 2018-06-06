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

#ifndef HIFLOW_EIGEN_VALUE_SLEPC_EIGENSOLVER_H_
#    define HIFLOW_EIGEN_VALUE_SLEPC_EIGENSOLVER_H_

#    include <cstdlib>
#    include "common/log.h"
#    include "eigen_solver.h"
#    include "mpi.h"
#    include "tools/mpi_tools.h"
#    include "common/sorted_array.h"
#    include "config.h"

#    ifdef WITH_PETSC
#        include "linear_algebra/petsc_la_descriptor.h"
#    endif

// TODO: Fix order dependency of next inclusion
#    include "common/smart_pointers.h"

namespace hiflow
{
    namespace la
    {

        /// Wrapper for performing matrix-free matrix-vector multiplication

        template <class LAD>
        struct Operator_wrapper
        {
            const typename LAD::MatrixType* op_;
            const typename LAD::MatrixType* opT_;
            const typename LAD::MatrixType* opi_;
            const typename LAD::MatrixType* opiT_;
            typename LAD::VectorType* vec_;

            int dof_offset_;
            int local_dofs_;
        };

        /// Forwarding PETSc helpers
        namespace petsc
        {
            class Mat_wrapper;
            class Vec_wrapper;
            class KSP_wrapper;
            class Val_wrapper;

            // Matrix type

            enum PmatType
            {
                MPI_SPARSE, // all matrices are explicitly given
                MPI_MATFREE, // all matrices are implicitly given 	
                MPI_MIXED // both explicit and implict matrices
            };
        }

        /// Forwarding SLEPc helpers
        namespace slepc
        {

            class EPS_wrapper;
            class SVD_wrapper;
            class ST_wrapper;

            struct MatrixStructure
            {
                /// Type of Matrix (MPI_SPARSE or MPI_MATFREE)
                petsc::PmatType mat_type;

                /// Column indices of nonzeros in diagonal part, obtained from SparsityStructure
                std::vector<int> diag_cols;
                /// Rows indices of nonzeros in diagonal part, obtained from SparsityStructure
                std::vector<int> diag_rows;

                /// Column indices of nonzeros in offdiagonal part, obtained from SparsityStructure
                std::vector<int> off_diag_cols;

                /// Row indices of nonzeros in offdiagonal part, obtained from SparsityStructure
                std::vector<int> off_diag_rows;

                /// Number of nonzeros per row in diagonal part
                std::vector<int> d_nnz;

                /// Number of nonzeros per row in offdiagonal part
                std::vector<int> o_nnz;

                /// Nonzero column indices for each row in diagonal block 
                std::vector< std::vector<int> > nz_col_ind_in_diag_row;

                /// Nonzero column indices for each row in offdiagonal block 
                std::vector< std::vector<int> > nz_col_ind_in_offdiag_row;

                /// Number of lcoal dofs
                int local_dofs;

                /// Number of nonzeros in diagonal part 
                int diag_nnz;

                /// Number of nonzeros in offdiagonal part
                int off_diag_nnz;

                /// Local dof offset 
                int dof_offset;

                /// Flag if struct is set up 
                bool is_set;
            };
        }

        /// @brief Base class for all eigenvalue / singular value solvers that make use of SLEPc

        template <class LAD>
        class SLEPcEigenSolver : public EigenSolver<LAD>
        {
          public:
            typedef typename LAD::MatrixType OperatorType;
            typedef typename LAD::VectorType VectorType;
            typedef typename LAD::DataType DataType;

            /// Constructor for use with EPS (EigenProblemSolver) object 
            /// @param[in] comm MPI communicator 
            /// @param[in] ptr_eps_wrapper Pointer to EPS wrapper object 
            SLEPcEigenSolver ( slepc::EPS_wrapper* ptr_eps_wrapper );

            /// Constructor for use with SVD (SingularValueDecomposition) object 
            /// @param[in] comm MPI communicator 
            /// @param[in] ptr_eps_wrapper Pointer to SVD wrapper object 
            SLEPcEigenSolver ( slepc::SVD_wrapper* ptr_svd_wrapper );

            /// Destructor

            virtual ~SLEPcEigenSolver ( )
            {
            };

            /// Setup operators for Eigenvalue / SingaluarValue problem
            /// @param[in] opA real part of left system matrix for eigenvalue problem A*x = k*B*x / matrix for SVD 
            /// @param[in] opAi imaginary part of left system matrix for eigenvalue problem A*x = k*B*x / matrix for SVD. Set to NULL, if no imaginary part is given 
            /// @param[in] cp Coupling object of considered matrix 
            /// @param[in] d_rows SparsityStrute.diagonal_rows 
            /// @param[in] d_cols SparsityStrute.diagonal_cols
            /// @param[in] o_rows SparsityStrute.off_diagonal_rows 
            /// @param[in] o_cols SparsityStrute.off_diagonal_cols
            void SetupOperatorA ( const OperatorType* opA, const OperatorType* opAi,
                                  const Couplings<DataType>& cp, const std::vector<int>& d_rows, const std::vector<int>& d_cols, const std::vector<int>& o_rows, const std::vector<int>& o_cols );

            /// Setup operator for Eigenvalue / SingaluarValue problem with matrix free implementation.
            /// Note: not every eigensolver needs the transposed operator. When using such a solver, use opAt = opA.
            /// @param[in] opA real part of left system matrix for eigenvalue problem A*x = k*B*x / matrix for SVD 
            /// @param[in] opAt real part of transposed left system matrix for eigenvalue problem A*x = k*B*x / matrix for SVD 
            /// @param[in] opAi (optional) imag part of left system matrix for eigenvalue problem A*x = k*B*x / matrix for SVD 
            /// @param[in] opAit (optional) imag part of transposed left system matrix for eigenvalue problem A*x = k*B*x / matrix for SVD 
            void SetupOperatorA ( const OperatorType* opA, const OperatorType* opAt, const OperatorType* opAi, const OperatorType* opAit, const Couplings<DataType>& cp );

            /// Setup operator for Eigenvalue / SingaluarValue problem
            /// @param[in] opB real part of right system matrix for eigenvalue problem A*x = k*B*x / matrix for SVD		
            /// @param[in] opBi (optional) imaginary part of right system matrix for eigenvalue problem A*x = k*B*x / matrix for SVD. Set to NULL, if no imaginary part is given 
            /// @param[in] cp Coupling object of considered matrix 
            /// @param[in] d_rows SparsityStrute.diagonal_rows 
            /// @param[in] d_cols SparsityStrute.diagonal_cols
            /// @param[in] o_rows SparsityStrute.off_diagonal_rows 
            /// @param[in] o_cols SparsityStrute.off_diagonal_cols
            void SetupOperatorB ( const OperatorType* opB, const OperatorType* opBi,
                                  const Couplings<DataType>& cp, const std::vector<int>& d_rows, const std::vector<int>& d_cols, const std::vector<int>& o_rows, const std::vector<int>& o_cols );

            /// Setup operator for Eigenvalue / SingaluarValue problem with matrix free implementation.
            /// Note: not every eigensolver needs the transposed operator. When using such a solver, use opAt = opA.
            /// @param[in] opB real part of right system matrix for eigenvalue problem A*x = k*B*x / matrix for SVD 
            /// @param[in] opBt real part of transposed right system matrix for eigenvalue problem A*x = k*B*x / matrix for SVD 
            /// @param[in] opBi (optional) imag part of right system matrix for eigenvalue problem A*x = k*B*x / matrix for SVD 
            /// @param[in] opBit (optional) imag part of transposed right system matrix for eigenvalue problem A*x = k*B*x / matrix for SVD 
            void SetupOperatorB ( const OperatorType* opB, const OperatorType* opBt, const OperatorType* opBi, const OperatorType* opBit, const Couplings<DataType>& cp );

            /// Call this functions after setting all operators defining the eigenvalue problem	
            virtual void PassOperatorsToSolver ( ) = 0;

            /// Initialize internal linear algebra structures used by solver
            /// @param[in] comm MPi communicator of considered matrix 

            virtual void Init ( const MPI_Comm& comm )
            {
            };

            /// Solves an Eigenvalue problem
            /// @return status if solver succeeded
            virtual EigenSolverState Solve ( ) = 0;

            /// Returns computed eigenvalues 
            /// @param[in] j Desired eigenvalue, sorting order corresponding to specified eigenvalue type
            /// @param[out] real Real part of j-th eigenvalue 
            /// @param[out] imag Imaginary part of j-th eigenvalue
            void GetEigenValue ( int j, DataType& real, DataType& imag );

            /// Return computed  eigenvectors
            /// @param[in] j Desired eigenvector, sorting order corresponding to specified eigenvalue type
            /// @param[out] real Real part of j-th eigenvector 
            /// @param[out] imag Imaginary part of j-th eigenvector
            void GetEigenVector ( int j, VectorType& v_real, VectorType& v_imag );

            /// Returns computed singular values 
            /// @param[in] j Desired singular value, sorting order corresponding to specified singular value type
            /// @param[out] sig j-th singular value 
            void GetSingularValue ( int j, DataType& sig );

            /// Return computed left and right singular vectors 
            /// @param[in] j Desired singular vector, sorting order corresponding to specified eigenvalue type
            /// @param[out] v_r j-th right singular vector 
            /// @param[out] v_l j-th left singular vector
            void GetSingularVector ( int j, VectorType& v_r, VectorType& v_l );

            /// Clear allocated data
            virtual void Clear ( );

          protected:
            /// Setup internal linear algebra structures used by solver
            /// @param[in] cp Coupling object of considered matrix 
            /// @param[in] d_rows SparsityStrute.diagonal_rows 
            /// @param[in] d_cols SparsityStrute.diagonal_cols
            /// @param[in] o_rows SparsityStrute.off_diagonal_rows 
            /// @param[in] o_cols SparsityStrute.off_diagonal_cols
            /// @param[out] m_struct MatrixStrucutre object
            void SetupMatrixStructure ( const Couplings<DataType>& cp, const std::vector<int>& d_rows, const std::vector<int>& d_cols, const std::vector<int>& o_rows, const std::vector<int>& o_cols,
                                        slepc::MatrixStructure& m_struct );

            /// Setup internal linear algebra structures used by solver
            /// @param[in] cp Coupling object of considered matrix 
            /// @param[out] m_struct MatrixStrucutre object
            void SetupMatrixStructure ( const Couplings<DataType>& cp, slepc::MatrixStructure& m_struct );

            /// Clear vector in matrix structure object 
            /// @param[in] m_struct
            void ClearMatrixStructure ( slepc::MatrixStructure& m_struct );

            /// Setup PETSC matrix object for operator in Eigenvalue / SingaluarValue problem
            /// @param[in] op (Hiflow) system matrix to be converted into PETSc matrix object
            /// @param[in] opi (optionla) (Hiflow) system matrix to be converted into PETSc matrix object, imag part
            /// @param[out] ptr_mat_wrapper PETSc matrix object
            void SetupPetscMat ( const OperatorType* op, const OperatorType* opi, slepc::MatrixStructure struct_real, slepc::MatrixStructure struct_imag, petsc::Mat_wrapper* ptr_mat_wrapper );

            /// Setup PETSC matrix object for operator in Eigenvalue / SingaluarValue problem
            /// @param[in] op (Hiflow) system matrix to be converted into PETSc matrix object
            /// @param[in] opT transposed (Hiflow) system matrix to be converted into PETSc matrix object
            /// @param[in] opi (optional) (Hiflow) system matrix to be converted into PETSc matrix object, imag part
            /// @param[in] opiT (optional) transposed (Hiflow) system matrix to be converted into PETSc matrix object, imag part
            /// @param[out] ptr_mat_wrapper PETSc matrix object
            void SetupPetscMatShell ( const OperatorType* op, const OperatorType* opT, const OperatorType* opi, const OperatorType* opiT, petsc::Mat_wrapper* ptr_mat_wrapper, int matrix );

            /// Init auxiliary vector for performing matrix vector product in matrx-free implementation
            void InitVec ( const Couplings<DataType>& cp );

            /// MPI rank;
            int rank_;

            /// Number of threads
            int nb_procs_;

            /// MPI communicator
            MPI_Comm comm_;

            /// SLEPc solver object eigenvalue problem solver object (EPS)
            hiflow::scoped_ptr<slepc::EPS_wrapper> ptr_eps_wrapper_;

            /// SLEPc solver object singular value decomposition object (SVD)
            hiflow::scoped_ptr<slepc::SVD_wrapper> ptr_svd_wrapper_;

            /// Pointer to PETSc matrix object for left system matrix
            petsc::Mat_wrapper* ptr_A_wrapper_;

            /// Pointer to PETSc matrix object for left system matrix
            petsc::Mat_wrapper* ptr_B_wrapper_;

            /// Wrapper object for performing matrix free matrix-vector multiplication with A 
            Operator_wrapper<LAD> opA_ptr_;

            /// Wrapper object for performing matrix free matrix-vector multiplication with B
            Operator_wrapper<LAD> opB_ptr_;

            /// Data Strucutures of operators 
            slepc::MatrixStructure struct_A_;
            slepc::MatrixStructure struct_Ai_;
            slepc::MatrixStructure struct_B_;
            slepc::MatrixStructure struct_Bi_;

            /// Real part of computed eigenvalues
            std::vector< petsc::Val_wrapper > eig_val_;

            /// Computed singular values 
            std::vector< DataType > sig_;

            /// complex computed eigenvectors / right singular vectors
            std::vector< petsc::Vec_wrapper* > ptr_v_wrapper_;

            /// Real part of computed eigenvectors / right singular vectors
            std::vector< petsc::Vec_wrapper* > ptr_vr_wrapper_;

            /// Imaginary part of computed eigenvectors / left singular vectors
            std::vector< petsc::Vec_wrapper* > ptr_vi_wrapper_;

            /// Auxiliary Hiflow vector for matrix-vector multiplication 
            VectorType aux_vec_;

            /// Flag indicating whether PETSC matrix A is initialized 
            bool mat_A_set_;

            /// Flag indicating whether PETSC matrix B is initialized 
            bool mat_B_set_;

        };

    } // namespace la
} // namespace hiflow

#endif  