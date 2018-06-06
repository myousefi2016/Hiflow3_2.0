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
#include "lpreconditioner_mc.h"
#include "lpreconditioner_mc_sgs.h"

#include <assert.h>

#include "lvector.h"
#include "lvector_cpu.h"

#include "lmatrix.h"
#include "lmatrix_csr_cpu.h"

namespace hiflow
{
    namespace la
    {

        template <typename ValueType>
        lPreconditioner_MultiColoring_SymmetricGaussSeidel<ValueType>::lPreconditioner_MultiColoring_SymmetricGaussSeidel ( )
        {
            this->precond_name_ = "Multi-Coloring / Symmetric Gauss-Seidel";
        }

        template <typename ValueType>
        lPreconditioner_MultiColoring_SymmetricGaussSeidel<ValueType>::~lPreconditioner_MultiColoring_SymmetricGaussSeidel ( )
        {
        }

        template <typename ValueType>
        void lPreconditioner_MultiColoring_SymmetricGaussSeidel<ValueType>::scale_D ( )
        {

            for ( int i = 0; i<this->ncolors_; ++i )
            {
                this->iDv_[i]->Scale ( -1.0 );
                this-> Dv_[i]->Scale ( -1.0 );

                // relaxation
                this->Dv_[i]->Scale ( 1.0 / ( this->omega_relax_ * ( 2.0 - this->omega_relax_ ) ) ); // for 1.0/(omega(2.0-omega))
                // note: the scaling of L/R is in the building procedure ( see extract_mat() )
            }

        }

        template <typename ValueType>
        void lPreconditioner_MultiColoring_SymmetricGaussSeidel<ValueType>::diagonalstep ( )
        {
            // Diag step
            // z:=-Dz
            for ( int i = 0; i<this->ncolors_; ++i )
                this->output_chunks_[i]->ElementWiseMult ( *this->Dv_[i] );
        }

        template <typename ValueType>
        void lPreconditioner_MultiColoring_SymmetricGaussSeidel<ValueType>::backwardstep ( )
        {

            int r_index = this->ncolors_ * ( this->ncolors_ - 1 ) / 2 - 1; // backward accessing

            for ( int i = this->ncolors_ - 1; i >= 0; --i )
            {

                // this->output_chunks_[i] += this->R_*this->output_chunks_[j]
                for ( int j = this->ncolors_ - 1; j > i; --j )
                {
                    this->R_[r_index]->VectorMultAdd ( *this->output_chunks_[j], this->output_chunks_[i] );
                    --r_index;
                }

                // this->output_chunks_[i] := -(D^-1)*this->output_chunks_[i];
                this->output_chunks_[i]->ElementWiseMult ( *this->iDv_[i] );
            }

        }

        template class lPreconditioner_MultiColoring_SymmetricGaussSeidel<double>;
        template class lPreconditioner_MultiColoring_SymmetricGaussSeidel<float>;

    } // namespace hiflow
} // namespace la
