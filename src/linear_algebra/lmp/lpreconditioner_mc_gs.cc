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
#include "lpreconditioner_mc_gs.h"

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
        lPreconditioner_MultiColoring_GaussSeidel<ValueType>::lPreconditioner_MultiColoring_GaussSeidel ( )
        {
            this->precond_name_ = "Multi-Coloring / Gauss-Seidel";
        }

        template <typename ValueType>
        lPreconditioner_MultiColoring_GaussSeidel<ValueType>::~lPreconditioner_MultiColoring_GaussSeidel ( )
        {
        }

        template <typename ValueType>
        void lPreconditioner_MultiColoring_GaussSeidel<ValueType>::Init ( )
        {
            assert ( false );
        }

        template <typename ValueType>
        void lPreconditioner_MultiColoring_GaussSeidel<ValueType>::build_aux_matrix ( )
        {
            assert ( this->Operator_ != NULL );
            this->aux_mat_analysis_cpu_ = new CPUsimple_CSR_lMatrix<ValueType>;
            this->aux_mat_analysis_cpu_->CloneFrom ( *this->Operator_ );
        }

        template <typename ValueType>
        void lPreconditioner_MultiColoring_GaussSeidel<ValueType>::allocate_LU ( )
        {
            this->LU_ = this->Operator_->CloneWithoutContent ( );
            this->LU_->CloneFrom ( *this->Operator_ );

            this->LU_cpu_ = this->aux_mat_analysis_cpu_->CloneWithoutContent ( );
            this->LU_cpu_->CloneFrom ( *this->Operator_ );

        }

        template <typename ValueType>
        void lPreconditioner_MultiColoring_GaussSeidel<ValueType>::factorize ( )
        {
            // do nothing ; matrix-free precond
        }

        template <typename ValueType>
        void lPreconditioner_MultiColoring_GaussSeidel<ValueType>::scale_D ( )
        {

            for ( int i = 0; i<this->ncolors_; ++i )
            {
                this->iDv_[i]->Scale ( -1.0 );
                this-> Dv_[i]->Scale ( -1.0 );

                // note: the scaling of L/R is in the building procedure ( see extract_mat() )
            }

        }

        template <typename ValueType>
        void lPreconditioner_MultiColoring_GaussSeidel<ValueType>::forwardstep ( )
        {

            int l_index = 0; // forward accessing

            for ( int i = 0; i<this->ncolors_; ++i )
            {

                // this->output_chunks_[i]:=-this->output_chunks_[i];
                this->output_chunks_[i]->Scale ( -1.0 );

                // this->output_chunks_[i] += this->L_*this->output_chunks_[j]
                for ( int j = 0; j < i; ++j )
                {
                    this->L_[l_index]->VectorMultAdd ( *this->output_chunks_[j], this->output_chunks_[i] );
                    ++l_index;
                }

                // this->output_chunks_[i] := -(D^-1)*this->output_chunks_[i];
                this->output_chunks_[i]->ElementWiseMult ( *this->iDv_[i] );
            }

        }

        template <typename ValueType>
        void lPreconditioner_MultiColoring_GaussSeidel<ValueType>::diagonalstep ( )
        {
            // Only scaling
            // z:=omega * z
            // note: the scaling of L/R is in the building procedure ( see extract_mat() )

            if ( this->omega_relax_ != 1.0 )
                for ( int i = 0; i<this->ncolors_; ++i )
                {
                    //      this->output_chunks_[i]->Scale(1.0/(this->omega_relax_));
                    this->output_chunks_[i]->Scale ( this->omega_relax_ );
                }

        }

        template <typename ValueType>
        void lPreconditioner_MultiColoring_GaussSeidel<ValueType>::backwardstep ( )
        {
            // do nothing
        }

        template class lPreconditioner_MultiColoring_GaussSeidel<double>;
        template class lPreconditioner_MultiColoring_GaussSeidel<float>;

    } // namespace hiflow
} // namespace la
