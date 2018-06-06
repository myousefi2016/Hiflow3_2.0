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

#ifndef HIFLOW_LINEARALGEBRA_PCE_VECTOR_
#    define HIFLOW_LINEARALGEBRA_PCE_VECTOR_

#    include "polynomial_chaos/pc_tensor.h"
#    include "hypre_vector.h"
#    include "coupled_vector.h"

using namespace hiflow::polynomialchaos;

namespace hiflow
{
    namespace la
    {

        template<class DataType>
        class PCEVector
        {
          public:

            // constructor
            PCEVector ( );

            // destructor
            ~PCEVector ( );

            // initialize
            void Init ( const PCTensor& pctensor, const HypreVector<DataType>& mean_vector );
            void Init ( const PCTensor& pctensor, const MPI_Comm& comm, const LaCouplings& cp );

            // access member of mode_vector_
            HypreVector<DataType>& Mode ( const int& mode );
            const HypreVector<DataType>& GetMode ( const int& mode ) const;

            // return number of modes
            int nb_mode ( ) const;

            // return total level
            int total_level ( ) const;

            // return pctensor
            PCTensor GetPCTensor ( ) const;

            // Update
            void Update ( );
            void Update ( const int& mode );

            // Clear
            void Clear ( );

            // Zeros
            void Zeros ( );

            // Axpy
            void Axpy ( const PCEVector<DataType>& vec, const DataType& alpha );
            void Axpy ( const PCEVector<DataType>& vec, const DataType& alpha, const int& l );
            void Axpy ( const int& mode, const PCEVector<DataType>& vec, const DataType& alpha );
            void Axpy ( const int& mode, const HypreVector<DataType>& vec, const DataType& alpha );

            // Dot
            DataType Dot ( const PCEVector<DataType>& vec ) const;
            DataType Dot ( const int& mode, const HypreVector<DataType>& vec ) const;

            // ScaleAdd
            void ScaleAdd ( const PCEVector<DataType>& vec, const DataType& alpha );
            void ScaleAdd ( const int& mode, const HypreVector<DataType>& vec, const DataType& alpha );

            // Scale
            void Scale ( const DataType& alpha );

            // size
            int size_local ( ) const;
            int size_local ( const int& mode ) const;

            int size_global ( ) const;
            int size_global ( const int& mode ) const;

            int size_local_ghost ( ) const;
            int size_local_ghost ( const int& mode ) const;

            // Norm2
            DataType Norm2 ( ) const;

            // CloneFrom
            void CloneFrom ( const PCEVector<DataType>& vec );
            void CloneFrom ( const PCEVector<DataType>& vec, const int& l );
            void CloneFrom ( const int& mode, const HypreVector<DataType>& vec );

            // CopyFrom
            void CopyFrom ( const PCEVector<DataType>& vec );
            void CopyFrom ( const PCEVector<DataType>& vec, const int & l );
            void CopyFrom ( const int& mode, const HypreVector<DataType>& vec );

            // CopyFromWithoutGhost
            void CopyFromWithoutGhost ( const PCEVector<DataType>& vec );
            void CopyFromWithoutGhost ( const PCEVector<DataType>& vec, const int& l );
            void CopyFromWithoutGhost ( const int& mode, const HypreVector<DataType>& vec );

            // CloneFromWithoutContent
            void CloneFromWithoutContent ( const PCEVector<DataType>& vec );
            void CloneFromWithoutContent ( const PCEVector<DataType>& vec, const int& l );
            void CloneFromWithoutContent ( const int& mode, const HypreVector<DataType>& vec );

            // Write and Read vector to/frome HDF5 file
            void WriteHDF5 ( const std::string& filename,
                             const std::string& groupname,
                             const std::string& datasetname );
            void ReadHDF5 ( const std::string& filename,
                            const std::string& groupname,
                            const std::string& datasetname );

          private:
            // set of mode vectors
            std::vector< HypreVector<DataType> > mode_vector_;

            // pc tensor
            PCTensor pctensor_;

            // number of modes
            int nmode_;

            // total level
            int totlevel_;

            // intialized flag
            bool initialized_;

        };

    }
}

#endif
