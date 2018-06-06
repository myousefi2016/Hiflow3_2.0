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

/// @author Chandramowli Subramanian, Martin Wlotzka

#ifndef HIFLOW_LINEARALGEBRA_COUPLINGS_H_
#    define HIFLOW_LINEARALGEBRA_COUPLINGS_H_

#    include "la_couplings.h"
#    include "dof/dof_fem_types.h"
#    include "mpi.h"

namespace hiflow
{

    namespace doffem
    {
        template<class DataType>
        class DofPartition;
    }

    namespace la
    {

        /// @brief Contains information for CoupledMatrix/CoupledVector.
        ///
        /// Contains border, ghost data and information about matrix structure.
        /// Uses @c hiflow::doffem::DofPartition internally.

        template<class DataType>
        class Couplings : public LaCouplings
        {
          public:

            typedef doffem::DofID DofID;

            /// Standard constructor, sets initialization status to @c false.
            /// @c Couplings must be initialized with @c Couplings::InitializeCouplings before use.
            Couplings ( );

            /// Destructor
            ~Couplings ( );

            /// Sets MPI communicator and dof partition.
            /// @param comm_hf HiFlow communicator
            /// @param dp dof partition
            void Init ( const MPI_Comm& comm_hf, const doffem::DofPartition<DataType>& dp );

            /// Computes couplings, i.e. needed border and ghost data. Sets initialization
            /// status to @c true if succesful.
            /// @param rows_offdiagonal row indices for diagonal block
            /// @param cols_offdiagonal column indices for diagonal block
            void InitializeCouplings ( const std::vector<int>& rows_offdiag,
                                       const std::vector<int>& cols_offdiag );

            /// @return Dof partition
            //virtual const doffem::DofPartition<DataType>& dof_partition() const
            //    { return *(this->dof_partition_); }

            const doffem::DofPartition<DataType>* dof_partition ( ) const
            {
                return this->dof_partition_;
            }

            Couplings& operator= ( const Couplings& cp )
            {
                if ( this == &cp )
                {
                    return *this;
                }

                // Clone MPI communicator
                if ( this->comm_ != MPI_COMM_NULL )
                {
                    MPI_Comm_free ( &( this->comm_ ) );
                }

                MPI_Comm_dup ( cp.comm ( ), &( this->comm_ ) );
                assert ( this->comm_ != MPI_COMM_NULL );

                // Get my rank
                MPI_Comm_rank ( this->comm_, &( this->my_rank_ ) );

                // Copy map global2offdiag
                this->global2offdiag_ = cp.global2offdiag ( );

                // Copy border_indices
                this->border_indices_ = cp.border_indices_vec ( );

                // Copy border offsets
                this->border_offsets_ = cp.border_offsets_vec ( );

                // Copy ghost offsets
                this->ghost_offsets_ = cp.ghost_offsets_vec ( );

                // Copy global offsets
                this->global_offsets_ = cp.global_offsets ( );

                // Copy offdiag2global
                this->offdiag2global_ = cp.offdiag2global ( );

                // Set initialization flag
                this->initialized_ = cp.initialized ( );

                // Copy pointer to DofPartition
                this->dof_partition_ = cp.dof_partition ( );

                return *this;
            }
          private:
            // no implementation of copy constructor or assignement operator
            Couplings ( const Couplings& );

            const doffem::DofPartition<DataType>* dof_partition_;
        };

    } // namespace la
} // namespace hiflow

#endif
