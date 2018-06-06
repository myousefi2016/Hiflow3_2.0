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

/// @author Martin Wlotzka, Mareike Schmidtobreick

#ifndef HIFLOW_LINEARALGEBRA_PP_VECTOR_H_
#    define HIFLOW_LINEARALGEBRA_PP_VECTOR_H_

#    include <iostream>
#    include "linear_algebra/lmp/la_global.h"
#    include "dof/dof_partition.h"

#    include "mpi.h"

namespace hiflow
{
    namespace la
    {

        /// @brief Post-processing vector.

        template<class LAD>
        class PpVector
        {
          public:

            typedef doffem::DofID DofID;
            typedef typename LAD::VectorType VectorType;
            typedef typename LAD::DataType DataType;

            /// standard constructor
            PpVector ( );

            /// destructor
            virtual ~PpVector ( );

            /// Sets MPI communicator and DofPartition
            void Init ( const MPI_Comm& comm, const doffem::DofPartition<DataType>& dp, const VectorType& vec );

            /// Retrieve DoF-IDs and offsets
            void InitStructure ( );

            /// Exchange values with other processes
            void UpdateValues ( );

            /// Get global DoF-IDs and values_
            void GetDofsAndValues ( std::vector<DofID>& id, std::vector<DataType>& val ) const;

            /// Get the values
            void GetValues ( std::vector<DataType>& val ) const;

            /// Get selected values
            void GetValues ( const int* global_dof_ids, const int num_values, DataType* out_values ) const;

            /// Get value at array-index (NOT DoF-ID!)
            DataType at ( const int i ) const;

            void Clear ( );

            void Print ( std::ostream& os ) const;

            inline int size ( )
            {
                return values_.size ( );
            }

          private:

            bool initialized_;

            int my_rank_;
            int nb_procs_;
            int mpi_info_;

            MPI_Comm comm_;
            const doffem::DofPartition<DataType>* dof_partition_;
            const VectorType* vector_;

            std::vector<DofID> dof_ids_;
            SortedArray<DofID> sorted_dof_ids_;
            std::vector<int> offsets_;
            std::vector<DataType> values_;

            std::vector<int> nb_send_;
            std::vector<int> nb_recv_;
            std::vector<DofID> send_dofs_;
            std::vector<int> send_offsets_;
            std::vector<DataType> send_values_;

        };

        template<class LAD>
        std::ostream& operator<< ( std::ostream& os, const PpVector<LAD>& vec )
        {
            vec.Print ( os );
            return os;
        }

    } // namespace la
} // namespace hiflow

#endif  // HIFLOW_LINEARALGEBRA_PP_VECTOR_H_
