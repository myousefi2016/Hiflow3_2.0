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

#include "lpreconditioner.h"
#include "lmatrix.h"
#include "lmatrix_csr_cpu.h"
#include "lmp_log.h"

#include <assert.h>

namespace hiflow
{
    namespace la
    {

        // Class lPreconditioner

        template <typename ValueType>
        lPreconditioner<ValueType>::lPreconditioner ( )
        {
            this->Operator_ = NULL;
            this->precond_name_ = "";
            this->reuse_ = true;
            this->modified_op_ = false;
            this->state_ = false;
        }

        template <typename ValueType>
        lPreconditioner<ValueType>::~lPreconditioner ( )
        {
        }

        template <typename ValueType>
        void lPreconditioner<ValueType>::SetupOperator ( const hiflow::la::lMatrix<ValueType> &op )
        {
            this->Operator_ = &op;
            this->SetModifiedOperator ( true );
        }

        template <typename ValueType>
        void lPreconditioner<ValueType>::Build ( )
        {
            this->SetState ( true );
            this->SetModifiedOperator ( false );
        }

        template <typename ValueType>
        void lPreconditioner<ValueType>::print ( std::ostream &out ) const
        {

            LOG_INFO ( "lPreconditioner", this->precond_name_ );
        }

        template <typename ValueType>
        void lPreconditioner<ValueType>::PermuteBack ( hiflow::la::lMatrix<ValueType> *op ) const
        {
            // do nothing by default
        }

        template <typename ValueType>
        void lPreconditioner<ValueType>::PermuteBack ( hiflow::la::lVector<ValueType> *vec ) const
        {
            // do nothing by default
        }

        template <typename ValueType>
        void lPreconditioner<ValueType>::Clear ( )
        {
            this->Operator_ = NULL;
            this->precond_name_ = "";
            this->reuse_ = false;
            this->modified_op_ = false;
            this->state_ = false;
        }

        // Class lPreconditioner_Jacobi

        template <typename ValueType>
        lPreconditioner_Jacobi<ValueType>::lPreconditioner_Jacobi ( )
        : lPreconditioner<ValueType>( )
        {
            this->precond_name_ = "Jacobi";
        }

        template <typename ValueType>
        lPreconditioner_Jacobi<ValueType>::~lPreconditioner_Jacobi ( )
        {
        }

        template <typename ValueType>
        void lPreconditioner_Jacobi<ValueType>::Init ( )
        {
            assert ( false );
        }

        template <typename ValueType>
        void lPreconditioner_Jacobi<ValueType>::Init ( const hiflow::la::lVector<ValueType> &output )
        {
            assert ( &output != NULL );

            this->inv_D_ = output.CloneWithoutContent ( );

        }

        template <typename ValueType>
        void lPreconditioner_Jacobi<ValueType>::Clear ( )
        {
            this->inv_D_->Clear ( );
            delete this->inv_D_;
            lPreconditioner<ValueType>::Clear ( );
        }

        template <typename ValueType>
        void lPreconditioner_Jacobi<ValueType>::Build ( )
        {
            assert ( this->Operator_ != NULL );
            assert ( this->Operator_->get_num_row ( ) == this->Operator_->get_num_col ( ) );

            CPU_CSR_lMatrix<ValueType> *mat_cpu = new CPUsimple_CSR_lMatrix<ValueType>;

            // Permute the Operator
            mat_cpu->CloneFrom ( *this->Operator_ );

            // Extract the inv-diagonal entries
            mat_cpu->extract_invdiagelements ( 0,
                                               mat_cpu->get_num_row ( ),
                                               this->inv_D_ );

            mat_cpu->Clear ( );

            delete mat_cpu;

            this->SetState ( true );
            this->SetModifiedOperator ( false );
        }

        template <typename ValueType>
        void lPreconditioner_Jacobi<ValueType>::ApplylPreconditioner ( const hiflow::la::lVector<ValueType> &input,
                                                                       hiflow::la::lVector<ValueType> *output )
        {
            if ( !this->GetState ( ) )
            {
                this->Build ( );
            }

            // CPU sequential
            //  this->Operator_->Pjacobi(input, output);

            // Parallel
            output->CopyFrom ( input );
            output->ElementWiseMult ( *this->inv_D_ );
        }

        // Class lPreconditioner_GaussSeidel

        template <typename ValueType>
        lPreconditioner_GaussSeidel<ValueType>::lPreconditioner_GaussSeidel ( )
        : lPreconditioner<ValueType>( )
        {
            this->precond_name_ = "Gauss-Seidel";
        }

        template <typename ValueType>
        lPreconditioner_GaussSeidel<ValueType>::~lPreconditioner_GaussSeidel ( )
        {
        }

        template <typename ValueType>
        void lPreconditioner_GaussSeidel<ValueType>::Init ( )
        {
            // do nothing - it's a matrix free preconditioner
        }

        template <typename ValueType>
        void lPreconditioner_GaussSeidel<ValueType>::ApplylPreconditioner ( const hiflow::la::lVector<ValueType> &input,
                                                                            hiflow::la::lVector<ValueType> *output )
        {
            assert ( this->Operator_ != NULL );
            if ( !this->GetState ( ) )
            {
                this->Build ( );
            }

            this->Operator_->Pgauss_seidel ( input, output );
        }

        // Class lPreconditioner_SymmetricGaussSeidel

        template <typename ValueType>
        lPreconditioner_SymmetricGaussSeidel<ValueType>::lPreconditioner_SymmetricGaussSeidel ( )
        : lPreconditioner<ValueType>( )
        {
            this->precond_name_ = "Symmetric Gauss-Seidel";
        }

        template <typename ValueType>
        lPreconditioner_SymmetricGaussSeidel<ValueType>::~lPreconditioner_SymmetricGaussSeidel ( )
        {
        }

        template <typename ValueType>
        void lPreconditioner_SymmetricGaussSeidel<ValueType>::Init ( )
        {
            // do nothing - it's a matrix free preconditioner
        }

        template <typename ValueType>
        void lPreconditioner_SymmetricGaussSeidel<ValueType>::ApplylPreconditioner ( const hiflow::la::lVector<ValueType> &input,
                                                                                     hiflow::la::lVector<ValueType> *output )
        {
            if ( !this->GetState ( ) )
            {
                this->Build ( );
            }

            this->Operator_->Psgauss_seidel ( input, output );
        }

        // Class lPreconditioner_BlocksSymmetricGaussSeidel

        template <typename ValueType>
        lPreconditioner_BlocksSymmetricGaussSeidel<ValueType>::lPreconditioner_BlocksSymmetricGaussSeidel ( )
        : lPreconditioner<ValueType>( )
        {
            this->precond_name_ = "Block-wise Symmetric Gauss-Seidel";
            this->num_blocks_ = 1;
        }

        template <typename ValueType>
        lPreconditioner_BlocksSymmetricGaussSeidel<ValueType>::~lPreconditioner_BlocksSymmetricGaussSeidel ( )
        {
        }

        template <typename ValueType>
        void lPreconditioner_BlocksSymmetricGaussSeidel<ValueType>::Init ( const int num_blocks )
        {
            this->num_blocks_ = num_blocks;
        }

        template <typename ValueType>
        void lPreconditioner_BlocksSymmetricGaussSeidel<ValueType>::Init ( void )
        {
            this->Init ( 1 );
        }

        template <typename ValueType>
        void lPreconditioner_BlocksSymmetricGaussSeidel<ValueType>::ApplylPreconditioner ( const hiflow::la::lVector<ValueType> &input,
                                                                                           hiflow::la::lVector<ValueType> *output )
        {
            if ( !this->GetState ( ) )
            {
                this->Build ( );
            }

            this->Operator_->BlocksPsgauss_seidel ( input, output, this->num_blocks_ );
        }

        // Class lPreconditioner_SOR

        template <typename ValueType>
        lPreconditioner_SOR<ValueType>::lPreconditioner_SOR ( )
        : lPreconditioner<ValueType>( )
        {
            this->precond_name_ = "SOR";
            this->relax_parameter_ = 1.0;
        }

        template <typename ValueType>
        lPreconditioner_SOR<ValueType>::~lPreconditioner_SOR ( )
        {
        }

        template <typename ValueType>
        void lPreconditioner_SOR<ValueType>::Init ( const ValueType omega )
        {
            this->relax_parameter_ = omega;
        }

        template <typename ValueType>
        void lPreconditioner_SOR<ValueType>::Init ( void )
        {
            this->Init ( 1.0 );
        }

        template <typename ValueType>
        void lPreconditioner_SOR<ValueType>::ApplylPreconditioner ( const hiflow::la::lVector<ValueType> &input,
                                                                    hiflow::la::lVector<ValueType> *output )
        {
            if ( !this->GetState ( ) )
            {
                this->Build ( );
            }

            this->Operator_->Psor ( this->relax_parameter_,
                                    input, output );
        }

        // Class lPreconditioner_SSOR

        template <typename ValueType>
        lPreconditioner_SSOR<ValueType>::lPreconditioner_SSOR ( )
        : lPreconditioner<ValueType>( )
        {
            this->precond_name_ = "Symmetric SOR";
            this->relax_parameter_ = 1.0;
        }

        template <typename ValueType>
        lPreconditioner_SSOR<ValueType>::~lPreconditioner_SSOR ( )
        {
        }

        template <typename ValueType>
        void lPreconditioner_SSOR<ValueType>::Init ( const ValueType omega )
        {
            this->relax_parameter_ = omega;
        }

        template <typename ValueType>
        void lPreconditioner_SSOR<ValueType>::Init ( void )
        {
            this->Init ( 1.0 );
        }

        template <typename ValueType>
        void lPreconditioner_SSOR<ValueType>::ApplylPreconditioner ( const hiflow::la::lVector<ValueType> &input,
                                                                     hiflow::la::lVector<ValueType> *output )
        {
            if ( !this->GetState ( ) )
            {
                this->Build ( );
            }

            this->Operator_->Pssor ( this->relax_parameter_,
                                     input, output );
        }

        // Class lPreconditioner_ILUp

        template <typename ValueType>
        lPreconditioner_ILUp<ValueType>::lPreconditioner_ILUp ( )
        : lPreconditioner<ValueType>( )
        {
            this->precond_name_ = "ILU(p)";
            this->ilu_p_ = 0;
            this->LU_ = NULL;
        }

        template <typename ValueType>
        lPreconditioner_ILUp<ValueType>::~lPreconditioner_ILUp ( )
        {
        }

        template <typename ValueType>
        void lPreconditioner_ILUp<ValueType>::Init ( int ilu_p )
        {
            this->ilu_p_ = ilu_p;
            assert ( this->ilu_p_ >= 0 );
        }

        template <typename ValueType>
        void lPreconditioner_ILUp<ValueType>::Init ( void )
        {
            this->Init ( 0 );
        }

        template <typename ValueType>
        void lPreconditioner_ILUp<ValueType>::Clear ( )
        {
            this->LU_->Clear ( );
            lPreconditioner<ValueType>::Clear ( );
        }

        template <typename ValueType>
        void lPreconditioner_ILUp<ValueType>::Build ( )
        {
            // make LU_ the same as Operator
            this->LU_ = this->Operator_->CloneWithoutContent ( );

            // make LU cpu matrix
            lMatrix<ValueType> *LU_cpu = new CPUsimple_CSR_lMatrix<ValueType>;

            LU_cpu->CloneFrom ( *this->Operator_ );

            // factorize cpu matrix
            if ( this->ilu_p_ == 0 )
            {
                LU_cpu->ilu0 ( );
            }
            else
            {
                LU_cpu->ilup ( this->ilu_p_ );
            }

            this->LU_->CloneFrom ( *LU_cpu );

            delete LU_cpu;
        }

        template <typename ValueType>
        void lPreconditioner_ILUp<ValueType>::ApplylPreconditioner ( const hiflow::la::lVector<ValueType> &input,
                                                                     hiflow::la::lVector<ValueType> *output )
        {
            if ( !this->GetState ( ) )
            {
                this->Build ( );
            }

            this->LU_->ilu_solve ( input,
                                   output );
        }

        template class lPreconditioner<double>;
        template class lPreconditioner<float>;

        template class lPreconditioner_Jacobi<double>;
        template class lPreconditioner_Jacobi<float>;

        template class lPreconditioner_GaussSeidel<double>;
        template class lPreconditioner_GaussSeidel<float>;

        template class lPreconditioner_SymmetricGaussSeidel<double>;
        template class lPreconditioner_SymmetricGaussSeidel<float>;

        template class lPreconditioner_BlocksSymmetricGaussSeidel<double>;
        template class lPreconditioner_BlocksSymmetricGaussSeidel<float>;

        template class lPreconditioner_SOR<double>;
        template class lPreconditioner_SOR<float>;

        template class lPreconditioner_SSOR<double>;
        template class lPreconditioner_SSOR<float>;

        template class lPreconditioner_ILUp<double>;
        template class lPreconditioner_ILUp<float>;

    } // namespace hiflow
} // namespace la
