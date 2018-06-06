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

/// @author Simon Gawlok

#ifndef HIFLOW_LINEARSOLVER_HYPRE_PRECONDITIONER_H_
#    define HIFLOW_LINEARSOLVER_HYPRE_PRECONDITIONER_H_

#    include <iostream>
#    include <cmath>

#    include "config.h"
#    include "common/log.h"
#    include "preconditioner.h"

#    ifdef WITH_HYPRE
#        include "_hypre_utilities.h"
#        include "HYPRE_krylov.h"
#        include "HYPRE.h"
#        include "HYPRE_parcsr_ls.h"
#    endif

namespace hiflow
{
    namespace la
    {

        /// @author Simon Gawlok
        /// @brief Base class for all Hypre preconditioners and linear solvers in HiFlow.

        template<class LAD>
        class HyprePreconditioner : virtual public Preconditioner<LAD>
        {
          public:
            typedef typename LAD::MatrixType OperatorType;
            typedef typename LAD::VectorType VectorType;
            typedef typename LAD::DataType DataType;

            HyprePreconditioner ( ) : comm_ ( MPI_COMM_NULL )
            {
                this->print_level_ = 0;
                if ( this->comm_ != MPI_COMM_NULL )
                {
                    MPI_Comm_free ( &this->comm_ );
                }
                this->SetModifiedParam ( false );
                this->SetInitialized ( false );
                this->SetState ( false );
                this->reuse_ = true;
                this->modified_op_ = false;
                this->use_solver_op_ = false;
                this->op_ = NULL;
                this->is_critical_hypre_solver_ = false;
            }

            virtual ~HyprePreconditioner ( )
            {
            }

#    ifdef WITH_HYPRE
            /// Get pointer to solve function of preconditioner
            virtual HYPRE_PtrToSolverFcn get_solve_function ( ) = 0;

            /// Get pointer to setup function of preconditioner
            virtual HYPRE_PtrToSolverFcn get_setup_function ( ) = 0;

            /// Get hypre preconditioner object
            virtual HYPRE_Solver& get_solver ( ) = 0;

            /// Get hypre preconditioner object

            HYPRE_Solver* GetSolver ( )
            {
                return this->solver_;
            }
#    endif

            virtual void Clear ( )
            {
                this->SetInitialized ( false );
                this->SetModifiedParam ( false );
                this->SetState ( false );
                this->print_level_ = 0;
                if ( this->comm_ != MPI_COMM_NULL )
                {
                    MPI_Comm_free ( &this->comm_ );
                }
                this->op_ = NULL;
                this->modified_op_ = false;
                this->reuse_ = true;
                this->use_solver_op_ = false;
                this->is_critical_hypre_solver_ = false;
            }

            /// Create Hypre solver object 

            virtual void Init ( )
            {
            };

            virtual void SetModifiedParam ( bool flag )
            {
                this->modified_param_ = flag;
                if ( flag )
                    this->SetState ( false );
            }

            virtual void SetInitialized ( bool flag )
            {
                this->initialized_ = flag;
                if ( !flag )
                    this->SetState ( false );
            }

          protected:
            /// MPI communicator
            MPI_Comm comm_;

#    ifdef WITH_HYPRE
            /// Hypre solver object
            HYPRE_Solver solver_;
#    endif   
            /// Flag if parameters have changed
            bool modified_param_;

            /// Flag if solver is initialized
            bool initialized_;

        };

    } // namespace la
} // namespace hiflow

#endif  // HIFLOW_LINEARSOLVER_HYPRE_PRECONDITIONER_H_
