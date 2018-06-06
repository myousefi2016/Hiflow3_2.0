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

#ifndef HIFLOW_POLYNOMIALCHAOS_PC_GALERKIN_VECTOR_H_
#    define HIFLOW_POLYNOMIALCHAOS_PC_GALERKIN_VECTOR_H_

/// \file pc_galerkin_vector.h
/// \brief Vector class for working with Galerkin projected Polynomial Chaos.
/// Each component of this PCGalerkinVector corresponds to a Finite-Element
/// vector related to some mode of the Polynomial Chaos expansion.
/// \author Michael Schick

#    include <vector>
#    include <cassert>
#    include <mpi.h>
#    include "linear_algebra/coupled_vector.h"
#    include "polynomial_chaos/pc_tensor.h"

namespace hiflow
{
    namespace polynomialchaos
    {

        template<class DataType>
        class PCGalerkinVector
        {
          public:
            typedef PCTensor GalerkinTensor;

            /// Default constructor
            PCGalerkinVector ( );
            /// Copy constructor
            PCGalerkinVector ( PCGalerkinVector<DataType> const& vec );
            /// Default destructor
            ~PCGalerkinVector ( );

            /// Clearing allocated memory
            void Clear ( );

            /// Returns pointer to mode with index n
            la::CoupledVector<DataType>* Mode ( int n ) const;
            /// Returns number of modes, which are owned by this MPI process
            inline int NModes ( ) const;
            /// Returns number of modes for the ghost layer
            inline int NOffModes ( ) const;

            /// Initialization of Vector by setting all modes
            void SetModes ( std::vector<la::CoupledVector<DataType>* > const& modes, int nmodes, int noff_modes = 0 );

            /// Axpy for specific mode index
            void ModeAxpy ( int mode, la::CoupledVector<DataType> const& vec, const DataType alpha );
            /// Set every mode to zero
            void Zeros ( );
            /// Axpy for two PCGalerkinVectors (all mode axpy)
            void Axpy ( const PCGalerkinVector<DataType>& vec, const DataType alpha );
            /// l^2 dot product between two PCGalerkinVectors
            DataType Dot ( const PCGalerkinVector<DataType>& vec ) const;
            /// Scale this with alpha and add vector vec
            void ScaleAdd ( const PCGalerkinVector<DataType>& vec, const DataType alpha );
            /// Scaling of the PCGalerkinVector (all modes)
            void Scale ( const DataType alpha );
            /// l^2 Norm of PCGalerkinVector
            DataType Norm2 ( ) const;

            /// Copy specific mode
            void CopyMode ( int mode, la::CoupledVector<DataType> const& vec );
            /// Clone complete PCGalerkinVector
            void CloneFrom ( const PCGalerkinVector<DataType>& vec );
            /// Copy complete PCGalerkinVector
            void CopyFrom ( const PCGalerkinVector<DataType>& vec );
            /// Clone structure of complete PCGalerkinVector
            void CloneFromWithoutContent ( const PCGalerkinVector<DataType>& vec );

            /// Local size of deterministic vector (Default mode index 0)

            inline int size_local ( ) const
            {
                return DSizeLocal ( );
            }
            /// Global size of deterministic vector (Default mode index 0)

            inline int size_global ( ) const
            {
                return DSizeLocal ( );
            }

            /// Update

            inline void Update ( ) const
            {
                this->UpdateCouplings ( );
            }
            /// Update Couplings of each mode

            inline void UpdateCouplings ( ) const
            {
                for ( int mode = 0; mode < nmodes_; ++mode )
                    this->modes_[mode]->Update ( );
            }
            /// Receive Ghost layer values (w.r.t. Finite-Element mesh)

            inline void ReceiveGhost ( )
            {
                for ( int mode = 0; mode < nmodes_; ++mode ) this->modes_[mode]->ReceiveGhost ( );
            }
            /// Send Ghost layer values (w.r.t. Finite-Element mesh)

            inline void SendBorder ( )
            {
                for ( int mode = 0; mode < nmodes_; ++mode ) this->modes_[mode]->SendBorder ( );
            }
            /// Wait until all MPI processors have received their mesh ghost layer

            inline void WaitForRecv ( )
            {
                for ( int mode = 0; mode < nmodes_; ++mode ) this->modes_[mode]->WaitForRecv ( );
            }
            /// Wait until all MPI processors have sent their mesh ghost layer

            inline void WaitForSend ( )
            {
                for ( int mode = 0; mode < nmodes_; ++mode ) this->modes_[mode]->WaitForSend ( );
            }
            /// Local size (Finite-Element dofs) of PC mode (Default mode index 0)
            inline int DSizeLocal ( ) const;
            /// Local size (Finite-Element dofs) of PC mode of specific mode index
            inline int DSizeLocal ( int mode ) const;

            /// Get values stored in deterministic vector of specific mode index
            inline void GetValues ( int mode, DataType* values ) const;
            /// Set values of deterministic vector of specific mode index
            inline void SetValues ( int mode, const DataType* values );

            /// Access to MPI communicator

            const MPI_Comm& comm ( ) const
            {
                return comm_;
            }
            /// Set MPI communicator
            void SetMPIComm ( const MPI_Comm& comm );

            /// MPI synchronization of ghost layer of Polynomial Chaos expansion
            void SynchronizeModes ( GalerkinTensor* pctensor );
            /// MPI synchronization of ghost layer of mesh ghost layer for all PC modes
            void SynchronizeModesGhost ( GalerkinTensor* pctensor );
            /// Change number of active modes owned by the processor

            void SetNLocalModes ( int nmodes )
            {
                nmodes_ = nmodes;
            }
            /// Coarsening of vector by restriction of polynomial degree of Chaos Polynomials
            void CreateCoarseVector ( GalerkinTensor* pctensor, PCGalerkinVector<DataType>& coarse_vec ) const;
            /// Refining vector by increasing polynomial degree of Chaos Polynomials
            void AssignToRefinedVector ( GalerkinTensor* pctensor, const PCGalerkinVector<DataType>& coarse_vec );

          protected:
            /// Basis initialization
            void Initialize ( );

            /// Number of owned mode, number of ghost layer modes, total number of modes and deterministic size
            int nmodes_, noff_modes_, nmodes_total_, size_local_deterministic_;
            /// MPI communicator
            MPI_Comm comm_;
            /// Container of externally allocated mode vectors
            std::vector<la::CoupledVector<DataType>* > modes_;
        };

        template<class DataType>
        int PCGalerkinVector<DataType>::NModes ( ) const
        {
            return nmodes_;
        }

        template<class DataType>
        int PCGalerkinVector<DataType>::NOffModes ( ) const
        {
            return noff_modes_;
        }

        template<class DataType>
        int PCGalerkinVector<DataType>::DSizeLocal ( ) const
        {
            return modes_[0]->size_local ( );
        }

        template<class DataType>
        int PCGalerkinVector<DataType>::DSizeLocal ( int mode ) const
        {
            return modes_[mode]->size_local ( );
        }

        template<class DataType>
        void PCGalerkinVector<DataType>::GetValues ( int mode, DataType* values ) const
        {
            modes_[mode]->GetValues ( values );
        }

        template<class DataType>
        void PCGalerkinVector<DataType>::SetValues ( int mode, const DataType* values )
        {
            modes_[mode]->SetValues ( values );
        }

    }
}

#endif
