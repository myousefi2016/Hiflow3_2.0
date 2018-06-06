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

#ifndef __LPRECONDITIONER_H
#    define __LPRECONDITIONER_H

#    include <iostream>
#    include <stdlib.h>

#    include "lvector.h"
#    include "lmatrix.h"

namespace hiflow
{
    namespace la
    {

        /// @brief Provides the base class to the local preconditioners
        /// @author Dimitar Lukarski

        template <typename ValueType>
        class lPreconditioner
        {
          public:

            lPreconditioner ( );
            virtual ~lPreconditioner ( );

            /// Initialize the preconditioner
            virtual void Init ( void ) = 0;

            /// Clear the preconditioner
            virtual void Clear ( void );

            /// Setup the matrix operator for the preconditioner
            /// @param op - set the matrix operator
            virtual void SetupOperator ( const hiflow::la::lMatrix<ValueType> &op );

            virtual void PermuteBack ( hiflow::la::lMatrix<ValueType> *op ) const;
            virtual void PermuteBack ( hiflow::la::lVector<ValueType> *vec ) const;

            /// Build (internally) the preconditioning matrix
            virtual void Build ( void );

            /// Apply the preconditioners (Solve Mz=r)
            /// @param - input vector r
            /// @return - output vector z
            virtual void ApplylPreconditioner ( const hiflow::la::lVector<ValueType> &input,
                                                hiflow::la::lVector<ValueType> *output ) = 0; // not const - MC algorithm!

            /// Set state of the preconditioner, i.e. whether it is ready to use or not.
            /// In case of reuse_ == false, this function always sets the state to false

            void SetState ( bool state )
            {
                this->state_ = state;
                if ( !this->reuse_ )
                    this->state_ = false;
            }

            /// Get State of the preconditioner

            bool GetState ( )
            {
                return this->state_;
            }

            /// Set flag whether preconditioner should be resued, in case it has not been changed
            /// Usually, this option should always be set to true. However, there might be MPI communicator problems when reusing
            /// too many BoomerAMG preconditioners at the same time

            void SetReuse ( bool flag )
            {
                this->reuse_ = flag;
                if ( !flag )
                    this->state_ = false;
            }

            bool GetReuse ( )
            {
                return this->reuse_;
            }

            /// Set status of operator

            void SetModifiedOperator ( bool flag )
            {
                this->modified_op_ = flag;
                if ( flag )
                    this->SetState ( false );
            }

            /// Print the type of preconditioner
            virtual void print ( std::ostream &out = std::cout ) const;

          protected:

            /// pointer to the operator matrix
            const hiflow::la::lMatrix<ValueType> *Operator_;

            std::string precond_name_;

            /// Flag if operator has changed
            bool modified_op_;

            /// Flag if preconditioner is set up. This flag is set to false, if either the operator or some parameters have changed
            bool state_;

            /// Flag if preconditioner should be reused, in case no changes have been made
            bool reuse_;

        };

        /// @brief Jacobi local Preconditioner
        /// @author Dimitar Lukarski

        template <typename ValueType>
        class lPreconditioner_Jacobi : public hiflow::la::lPreconditioner<ValueType>
        {
          public:

            lPreconditioner_Jacobi ( );
            virtual ~lPreconditioner_Jacobi ( );

            virtual void Init ( void );
            virtual void Init ( const hiflow::la::lVector<ValueType> &output );
            virtual void Clear ( void );

            virtual void Build ( void );

            virtual void ApplylPreconditioner ( const hiflow::la::lVector<ValueType> &input,
                                                hiflow::la::lVector<ValueType> *output );

          protected:

            hiflow::la::lVector<ValueType> *inv_D_;

        };

        /// @brief Gauss Seidel local Preconditioner
        /// @author Dimitar Lukarski

        template <typename ValueType>
        class lPreconditioner_GaussSeidel : public hiflow::la::lPreconditioner<ValueType>
        {
          public:

            lPreconditioner_GaussSeidel ( );
            virtual ~lPreconditioner_GaussSeidel ( );

            virtual void Init ( void );

            virtual void ApplylPreconditioner ( const hiflow::la::lVector<ValueType> &input,
                                                hiflow::la::lVector<ValueType> *output );

        };

        /// @brief Symmetric Gauss Seidel local Preconditioner
        /// @author Dimitar Lukarski

        template <typename ValueType>
        class lPreconditioner_SymmetricGaussSeidel : public hiflow::la::lPreconditioner<ValueType>
        {
          public:

            lPreconditioner_SymmetricGaussSeidel ( );
            virtual ~lPreconditioner_SymmetricGaussSeidel ( );

            virtual void Init ( void );

            virtual void ApplylPreconditioner ( const hiflow::la::lVector<ValueType> &input,
                                                hiflow::la::lVector<ValueType> *output );

        };

        /// @brief Block Symmetric Gauss-Seidel local Preconditioner
        /// - each of the block is executed in parallel according to the platform implementation
        /// @author Dimitar Lukarski

        template <typename ValueType>
        class lPreconditioner_BlocksSymmetricGaussSeidel : public hiflow::la::lPreconditioner<ValueType>
        {
          public:

            lPreconditioner_BlocksSymmetricGaussSeidel ( );
            virtual ~lPreconditioner_BlocksSymmetricGaussSeidel ( );

            virtual void Init ( void );
            virtual void Init ( const int num_blocks );

            virtual void ApplylPreconditioner ( const hiflow::la::lVector<ValueType> &input,
                                                hiflow::la::lVector<ValueType> *output );

          private:

            int num_blocks_;

        };

        /// @brief SOR local Preconditioner
        /// @author Dimitar Lukarski

        template <typename ValueType>
        class lPreconditioner_SOR : public hiflow::la::lPreconditioner<ValueType>
        {
          public:

            lPreconditioner_SOR ( );
            virtual ~lPreconditioner_SOR ( );

            virtual void Init ( void );
            virtual void Init ( const ValueType omega );

            virtual void ApplylPreconditioner ( const hiflow::la::lVector<ValueType> &input,
                                                hiflow::la::lVector<ValueType> *output );

          protected:

            ValueType relax_parameter_;

        };

        /// @brief Symmetric SOR local Preconditioner
        /// @author Dimitar Lukarski

        template <typename ValueType>
        class lPreconditioner_SSOR : public hiflow::la::lPreconditioner<ValueType>
        {
          public:

            lPreconditioner_SSOR ( );
            virtual ~lPreconditioner_SSOR ( );

            virtual void Init ( void );
            virtual void Init ( const ValueType omega );

            virtual void ApplylPreconditioner ( const hiflow::la::lVector<ValueType> &input,
                                                hiflow::la::lVector<ValueType> *output );

          protected:

            ValueType relax_parameter_;

        };

        /// @brief Symmetric ILU(p) local Preconditioner
        /// @author Dimitar Lukarski

        template <typename ValueType>
        class lPreconditioner_ILUp : public hiflow::la::lPreconditioner<ValueType>
        {
          public:

            lPreconditioner_ILUp ( );
            virtual ~lPreconditioner_ILUp ( );

            virtual void Init ( void );
            virtual void Init ( int ilu_p );
            virtual void Clear ( void );

            virtual void Build ( void );

            virtual void ApplylPreconditioner ( const hiflow::la::lVector<ValueType> &input,
                                                hiflow::la::lVector<ValueType> *output );

          protected:

            int ilu_p_;
            hiflow::la::lMatrix<ValueType> *LU_;

        };

    } // namespace la
} // namespace hiflow

#endif
