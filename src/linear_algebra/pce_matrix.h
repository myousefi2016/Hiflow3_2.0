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

#ifndef HIFLOW_LINEARALGEBRA_PCE_MATRIX_H_
#    define HIFLOW_LINEARALGEBRA_PCE_MATRIX_H_

#    include "polynomial_chaos/pc_tensor.h"
#    include "hypre_matrix.h"
#    include "coupled_matrix.h"
#    include "pce_vector.h"

using namespace hiflow::polynomialchaos;

namespace hiflow
{
    namespace la
    {

        template<class DataType>
        class PCEMatrix
        {
          public:

            // constructor
            PCEMatrix ( );

            // destructor
            ~PCEMatrix ( );

            // Inititialize
            void Init ( PCTensor& pctensor, const MPI_Comm& comm, const LaCouplings& cp );

            // define operator to access member of basis_matrix_
            HypreMatrix<DataType>& BasisMode ( const int& i );
            const HypreMatrix<DataType>& GetBasisMode ( const int& i ) const;

            // number of basis matrices
            int nb_basis ( ) const;

            // Zeros
            void Zeros ( );
            void Zeros ( const int& i );

            // Clear
            void Clear ( );

            // VectorMult
            void VectorMult ( PCEVector<DataType>& in, PCEVector<DataType> *out ) const;
            void VectorMult ( PCEVector<DataType>& in, PCEVector<DataType> *out, const int& l ) const;
            void VectorMult ( const int& i, HypreVector<DataType>& in, HypreVector<DataType> *out ) const;

          private:

            // a set of basis matrices
            std::vector< HypreMatrix<DataType> > basis_matrix_;

            // pc tensor
            PCTensor pctensor_;

            // size of basis
            int nbasis_;
        };

    }
}

#endif // HIFLOW_LINEARALGEBRA_PCE_MATRIX_H_
