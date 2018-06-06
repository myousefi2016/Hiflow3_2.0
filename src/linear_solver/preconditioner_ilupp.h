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

/// \author Chandramowli Subramanian, Hendryk Bockelmann

#ifndef HIFLOW_LINEARSOLVER_PRECONDITIONER_ILUPP_H_
#    define HIFLOW_LINEARSOLVER_PRECONDITIONER_ILUPP_H_

#    include <vector>

#    include "config.h"
#    include "linear_solver/preconditioner_bjacobi.h"

#    ifdef WITH_ILUPP
#        include "iluplusplus_interface.h"
#    endif

namespace hiflow
{
    namespace la
    {

        /// \brief ILU++ preconditioner interface
        ///

        template<class LAD>
        class PreconditionerIlupp : public PreconditionerBlockJacobi<LAD>
        {
          public:
            typedef typename LAD::MatrixType OperatorType;
            typedef typename LAD::VectorType VectorType;
            typedef typename LAD::DataType DataType;

            /// standard constructor
            PreconditionerIlupp ( );
            /// destructor
            virtual ~PreconditionerIlupp ( );

            /// Inits parameters for ILU++ preconditioner.
            /// \param prepro_type type of preprocessing
            /// \param precond_no number of preconditioner
            /// \param max_levels maximum number of multilevels
            /// \param mem_factor see ILU++ manual
            /// \param threshold see ILU++ manual
            /// \param min_pivot see ILU++ manual
            void InitParameter ( int prepro_type, int precond_no,
                                 int max_levels, double mem_factor,
                                 double threshold, double min_pivot );

            /// Computes the incomplete LU factorization.
            void Factorize ( );

            /// Build the preconditioner
            void Build ( );

            /// Applies the ILU++ preconditioner.
            /// @param b right hand side vector
            /// @param x solution vector
            /// @return status if preconditioning succeeded
            LinearSolverState ApplyPreconditioner ( const VectorType& b, VectorType* x );

            /// Clears allocated data.
            void Clear ( );

          protected:
            /// Creates local matrix in CSR format.
            void CreateLocalMatrix ( );

            // local matrix in CSR format
            std::vector<int> ia_;
            std::vector<int> ja_;
            std::vector<DataType> val_;

#    ifdef WITH_ILUPP
            iluplusplus::iluplusplus_precond_parameter ilupp_param_;
            iluplusplus::multilevel_preconditioner ilupp_precond_;
            iluplusplus::preprocessing_sequence ilupp_preproc_;
#    endif

        };

    } // namespace la
} // namespace hiflow

#endif  // HIFLOW_LINEARSOLVER_PRECONDITIONER_ILUPP_H_
