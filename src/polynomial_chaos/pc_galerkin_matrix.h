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

#ifndef HIFLOW_POLYNOMIALCHAOS_PC_GALERKIN_MATRIX_H_
#    define HIFLOW_POLYNOMIALCHAOS_PC_GALERKIN_MATRIX_H_

/// \file pc_galerkin_matrix.h
/// \brief Matrix class for working with Galerkin projected Polynomial Chaos.
/// using PCGalerkinVector class. The modes of the matrix correspond to the modes
/// of a Polynomial Chaos expansion of the differential operator modelling the PDE.
/// \author Michael Schick

#    include "linear_algebra/coupled_matrix.h"
#    include "polynomial_chaos/pc_tensor.h"
#    include "polynomial_chaos/pc_galerkin_vector.h"
#    include <vector>
#    include <mpi.h>

namespace hiflow
{
    namespace polynomialchaos
    {

        template<class DataType>
        class PCGalerkinMatrix
        {
          public:
            typedef PCTensor GalerkinTensor;

            /// Default constructor

            PCGalerkinMatrix ( )
            {
            }
            /// Constructor using a dummy deterministic vector for internal allocation
            PCGalerkinMatrix ( la::CoupledVector<DataType> const& dummy );
            /// Destructor
            ~PCGalerkinMatrix ( );

            /// For use with OpenMP: Set number of threads
            void SetNumThreads ( int numthreads );
            /// Set matrix mode entries with externally allocated memory
            void SetMatrix ( std::vector<la::CoupledMatrix<DataType>*> modes );
            /// Set matrix mode entires with internally allocated memory
            void SetMatrixAllocate ( std::vector<la::CoupledMatrix<DataType>*> modes );
            /// Set Polynomial Chaos Tensor
            void SetTensor ( GalerkinTensor* pctensor );
            /// Matrix-Vector multiplication
            void VectorMult ( PCGalerkinVector<DataType>& in,
                              PCGalerkinVector<DataType>* out ) const;
            /// Set every matrix mode to zero
            void Zeros ( );

            /// Clone complete matrix
            void CloneFrom ( const PCGalerkinMatrix<DataType>& mat );
            /// Copy complete matrix
            void CopyFrom ( const PCGalerkinMatrix<DataType>& mat );

            /// Returns a pointer to a specific mode of the Galerkin matrix

            la::CoupledMatrix<DataType>* Mode ( int mode ) const
            {
                return modes_[mode];
            }
            /// Returns pointer to a container of the matrix modes
            std::vector<la::CoupledMatrix<DataType>*> const* GetModes ( ) const;
            /// Returns pointer to Polynomial Chaos Tensor
            GalerkinTensor* GetTensor ( ) const;

            /// Number of lcoal rows of default mode index 0
            inline int nrows_local ( ) const;
            /// Number of local nonzero entries of default mode index 0
            inline int nnz_local_diag ( ) const;
            /// Extract CSR structure of default mode index 0
            inline void ExtractDiagonalCSR ( int* ia, int* ja, DataType* val ) const;

            /// Set level (total polynomial degree) for multilevel approach

            inline void SetLevel ( int level )
            {
                level_ = level;
                pctensor_->SetLevel ( level );
            }
            /// Get current level (total polynomial degree) for multilevel approach

            inline int GetCurrentLevel ( ) const
            {
                return level_;
            }

          protected:
            /// Number of modes for the matrix
            int nmodes_;
            /// Level = total polynomial degree of Chaos Polynomials
            int level_;
            /// Check variable if modes are stored externally or internally
            bool modes_allocated_;
            /// Container for matrix modes
            std::vector<la::CoupledMatrix<DataType>*> modes_;
            /// Polynomial Chaos Tensor
            GalerkinTensor* pctensor_;
            /// Number of OpenMP threads
            int numthreads_;
            /// Container for dummy variables for OpenMP parallelization
            std::vector<la::CoupledVector<DataType>*> tmp_vec_;
        };

        template<class DataType>
        int PCGalerkinMatrix<DataType>::nrows_local ( ) const
        {
            return modes_[0]->nrows_local ( );
        }

        template<class DataType>
        int PCGalerkinMatrix<DataType>::nnz_local_diag ( ) const
        {
            return modes_[0]->nnz_local_diag ( );
        }

        template<class DataType>
        void PCGalerkinMatrix<DataType>::ExtractDiagonalCSR ( int* ia, int* ja, DataType* val ) const
        {
            modes_[0]->ExtractDiagonalCSR ( ia, ja, val );
        }

    }
}

#endif
