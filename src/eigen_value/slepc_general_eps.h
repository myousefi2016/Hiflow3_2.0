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

#ifndef HIFLOW_EIGEN_VALUE_SLEPC_GENERAL_EPS_H_
#    define HIFLOW_EIGEN_VALUE_SLEPC_GENERAL_EPS_H_

#    include <mpi.h>
#    include "slepc_eigen_solver.h"

namespace hiflow
{
    namespace la
    {
        namespace slepc
        {

            enum EpsType
            {
                POWER, // Power method    
                SUBSPACE, // Subspace iteration
                ARNOLDI, // Arnoldi
                LANCZOS, // Lanczos
                KRYLOVSCHUR, // Krylov-Schur
                GD, // Generalized Davidson
                JD, // Jacobi-Davidson
                RQCG, // Rayleigh quotient CG
                ARPACK // Arpack wrapper
            };

            enum EpsProblemType
            {
                HEP, // Hermitian
                NHEP, // Non-Hermitian
                GHEP, // Generalized Hermitian
                GHIEP, // Generalized Hermitian indefinite
                GNHEP, // Generalized Non-Hermitian
                PGNHEP // GNHEP with positive (semi-) definite B  
            };

            enum EpsWhich
            {
                L_MAG, // Largest magnitude
                S_MAG, // Smallest magnitude
                L_REAL, // Largest real part     
                S_REAL, // Smallest real part 
                L_IMAG, // Largest imaginary part 
                S_IMAG, // Smallest imaginary part    
                ALL, // All eigenvalues in interval 
                TARGET // Eigenvalues closest to target value
            };

            enum EpsConv
            {
                ABS, // absolute value |r|
                REL, // relative to eigenvalue |r| / |nu|
                NORM // relative to matrix norm |r| / (|A| + |nu| |B|)
            };

            enum EpsMeasure
            {
                MAG, // complex difference 
                REAL, // difference in real part
                IMAG // difference in imaginary part    
            };

            enum StType
            {
                SHIFT, // B^{-1}A - sig I    
                SINVERT, // (A - sig B)^{-1} B
                CAYLEY, // (A - sig B)^{-1} (A + nu B) 
                PRECOND // K^{-1} ~ (A-sigB)^{-1}    
            };

            enum KSPType
            {
                GMRES,
                FGMRES,
                CG,
                PREONLY
            };

            enum PcType
            {
                LU,
                ILU,
                NONE,
                JACOBI,
                BJACOBI
            };

            enum SolverPackage
            {
                MUMPS,
                UMFPACK,
                NOPACKAGE
            };

            /// Forwarding Slepc EPS class
            class EPS_wrapper;
            class SVD_wrapper;
            class ST_wrapper;
        } // namespace slepc

        namespace petsc
        {
            class KSP_wrapper;
            class PC_wrapper;
        }

        /// @brief Wrapper class for general EPS implementations of SLEPc

        template <class LAD>
        class SLEPcGeneralEPS : public SLEPcEigenSolver<LAD>
        {
          public:
            typedef typename LAD::MatrixType OperatorType;
            typedef typename LAD::VectorType VectorType;
            typedef typename LAD::DataType DataType;

            /// Constructor
            /// Possible choices for EPS type 
            /// POWER      :  Power method    
            /// SUBSPACE   :  Subspace iteration
            /// ARNOLDI    :  Arnoldi
            ///    LANCZOS    :  Lanczos
            /// KRYLOVSCHUR: Krylov-Schur
            /// GD         : Generalized Davidson
            /// JD         : Jacobi-Davidson
            /// RQCG       : Rayleigh quotient CG
            /// ARPACK     : Arpack wrapper (Implictly restarted Arnoldi)
            /// Possible choices for mat_type: MPI_SPARSE 
            /// @param[in] eps_type Type of eigenvalue solver EPS 
            /// @param[in] mat_type Type of considered matrix
            SLEPcGeneralEPS ( slepc::EpsType eps_type );

            /// Destructor
            ~SLEPcGeneralEPS ( );

            /// Call this functions after setting all operators defining the eigenvalue problem    
            virtual void PassOperatorsToSolver ( );

            /// Initialize internal linear algebra structures used by solver
            /// @param[in] comm MPi communicator of considered matrix 
            virtual void Init ( const MPI_Comm& comm );

            /// Set dimensions, tolerances and convergenve criterion
            /// SLEPc documentation: " The parameters ncv and mpd are intimately related, so that the user is advised to set one of them at most. 
            /// Normal usage is that (a) in cases where nev is small, the user sets ncv (a reasonable default is 2*nev); and (b) in cases where nev is large, the user sets mpd.
            /// The value of ncv should always be between nev and (nev+mpd), typically ncv=nev+mpd. If nev is not too large, mpd=nev is a reasonable choice, otherwise a smaller value should be used." 
            /// Possible choices for conv_type and error_type:
            /// ABS : absolute value |r|
            /// REL : relative to eigenvalue |r| / |nu|
            /// NORM: relative to matrix norm |r| / (|A| + |nu| |B|)
            /// @param[in] nev          number of eigenvalues to be computed
            /// @param[in] ncv        largest dimension of working subspace, if set to 0: use SLEPc default 
            /// @param[in] mpd        maximum projected dimension, if set to 0: use SLEPc default 
            /// @param[in] maxits     maximum number of iterations
            /// @param[in] tol        residual tolerance 
            /// @param[in] conv_type  type of applied convergence criterion  
            /// @param[in] error_type type of error estimator 
            /// @param[in] true_res_norm use true (explicitly) calculated residual norm, instead of update formula, default: true    
            void InitControl ( int nev, int ncd, int mpd, int maxits, double tol, slepc::EpsConv conv_type, slepc::EpsConv error_type, bool true_res_norm );

            /// Set type of problem
            /// Possible choices for type:
            ///    HEP   : Hermitian
            /// NHEP  : Non-Hermitian
            /// GHEP  : Generalized Hermitian
            /// GHIEP : Generalized Hermitian indefinite
            /// GNHEP :    Generalized Non-Hermitian
            /// PGNHEP: GNHEP with positive (semi-) definite B  
            /// @param[in] type Type of problem
            void SetProblemType ( slepc::EpsProblemType type );

            /// Compute single eigenvalue, specified by type 
            /// Possible choices for type:  
            /// L_MAG :    Largest magnitude
            /// S_MAG :    Smallest magnitude
            /// L_REAL: Largest real part     
            /// S_REAL: Smallest real part 
            /// L_IMAG: Largest imaginary part 
            /// S_IMAG: Smallest imaginary part    
            /// @param[in] Type of extremal eigenvalue
            void ComputeExtremalEigenValue ( slepc::EpsWhich type );

            /// Compute eigenvalue closest to target
            /// Possible choices for type:
            ///    MAG : complex difference 
            /// REAL: difference in real part
            /// IMAG: difference in imaginary part    
            /// @param[in] target Compute eigenvalues closest to this point 
            /// @param[in] type Type of distance measure
            void ComputeTargetEigenValue ( DataType target, slepc::EpsMeasure type );

            /// Compute all eigenvalues in given interval  
            /// @param[in] a Lower bound
            /// @param[in] b Upper bound
            void ComputeIntervalEigenValue ( DataType a, DataType b );

            /// Set initial subspace @TODO 
            void SetInitialSpace ( std::vector< VectorType > is );

            /// Check whether chosen solver is applicable for defined problem 
            int CheckSolver ( );

            /// Solves an Eigenvalue problem 
            /// @return status if solver succeeded
            EigenSolverState Solve ( );

            /// Set shift of spectral transformation 
            void SetSTShift ( double shift );

            /// Set type of spectral transformation 
            /// Possible choices for type: (sig = shift, nu=estimate for eigenvalue)
            ///    SHIFT  : B^{-1}A - sig I    
            /// SINVERT: (A - sig B)^{-1} B
            /// CAYLEY : (A - sig B)^{-1} (A + nu B) 
            /// PRECOND: K^{-1} ~ (A-sigB)^{-1}    
            /// @param[in] type Type of spectral transformation
            void SetSTType ( slepc::StType type );

            /// Extract KSP solver object of spectral transformation
            petsc::KSP_wrapper* GetSTKSP ( );

            /// Extract PC preconditioner object of spectral transformation
            petsc::PC_wrapper* GetSTPC ( );

            /// Set Solver type of KSP in ST object 
            void SetSTKSP ( slepc::KSPType, int max_iter, int max_size, double abs_tol, double rel_tol );

            /// Set Preconditioner type in KSP object in ST object
            void SetSTPC ( slepc::PcType, slepc::SolverPackage );

            /// Setup Solver object of ST, needs to be called after parameters of solver are set
            void SetupSTSolver ( );

            /// @TODO add more functionalities to deal with ST object

            /// Clear allocated data 
            void Clear ( );

          private:
            /// SLEPc spectral transformation object (ST)
            hiflow::scoped_ptr<slepc::ST_wrapper> ptr_st_wrapper_;

            /// KSP solver object in ST
            petsc::KSP_wrapper* ptr_st_ksp_wrapper_;

            /// PC preconditioner object in ST 
            petsc::PC_wrapper* ptr_st_pc_wrapper_;

            /// Operator object in ST 
            petsc::Mat_wrapper* ptr_st_mat_wrapper_;

            /// Solver type
            slepc::EpsType eps_type_;

            /// Problem type 
            slepc::EpsProblemType problem_type_;

            /// Eigenvalue Type
            slepc::EpsWhich eval_type_;

            /// Error type 
            slepc::EpsConv error_type_;

            /// Convergence Type 
            slepc::EpsConv conv_type_;

        };

    } // namespace la
} // namespace hiflow

#endif  
