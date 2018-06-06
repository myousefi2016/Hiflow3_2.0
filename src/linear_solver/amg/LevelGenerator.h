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

/// @author Bernd Doser, HITS gGmbH
/// @date 2015-09-29

#ifndef SRC_LINEAR_SOLVER_AMG_LEVELGENERATOR_H_
#    define SRC_LINEAR_SOLVER_AMG_LEVELGENERATOR_H_

#    include "common/smart_pointers.h"

namespace hiflow
{
    namespace AMG
    {

        /**
         *  /brief Generation functor of a coarsening level.
         *
         *  To get a level hierarchy the LevelGenerator have to be called recursively,
         *  using the resulting coarse-matrix as input for the next step.
         *
         *  The LevelGeneratro is restricted on negatively coupled M-Matrices.
         */
        template <
        class MatrixArg,
        template <class> class CoarseningArg,
        template <class> class GridTransferArg,
        template <class> class CoarseMatrixConstructionArg
        >
        class LevelGenerator
        {
          public:

            // Due to the fact that the local and coupled matrices are not copyable or movable,
            // non-const matrices must be handled as SharedPtr, whereas const matrices will be
            // handled as const reference to avoid SharedPtr<const MatrixType>.
            typedef MatrixArg MatrixType;
            typedef hiflow::shared_ptr<MatrixType> PtrMatrixType;
            typedef LevelGenerator<MatrixType, CoarseningArg, GridTransferArg, CoarseMatrixConstructionArg> Self;

            typedef CoarseningArg<Self> CoarseningType;
            typedef GridTransferArg<Self> GridTransferType;
            typedef CoarseMatrixConstructionArg<Self> CoarseMatrixConstructionType;

            // ResultType must be defined before settings typedefs,
            // because ResultType is used within the methods.

            struct ResultType
            {
                PtrMatrixType ptr_interpolation_matrix;
                PtrMatrixType ptr_restriction_matrix;
                PtrMatrixType ptr_coarse_matrix;

                hiflow::shared_ptr<typename CoarseningType::Output> ptr_coarsening_output;
                hiflow::shared_ptr<typename GridTransferType::Output> ptr_grid_transfer_output;
                hiflow::shared_ptr<typename CoarseMatrixConstructionType::Output> ptr_coarse_matrix_construction_output;
            };

            typedef typename CoarseningType::Settings CoarseningSettings;
            typedef typename GridTransferType::Settings GridTransferSettings;
            typedef typename CoarseMatrixConstructionType::Settings CoarseMatrixConstructionSettings;

            LevelGenerator (
                             CoarseningSettings const& coarsening_settings = CoarseningSettings ( ),
                             GridTransferSettings const& grid_transfer_settings = GridTransferSettings ( ),
                             CoarseMatrixConstructionSettings const& coarse_matrix_construction_settings = CoarseMatrixConstructionSettings ( )
                             ) :
            coarsening_settings_ ( coarsening_settings ),
            grid_transfer_settings_ ( grid_transfer_settings ),
            coarse_matrix_construction_settings_ ( coarse_matrix_construction_settings )
            {
            }

            ResultType operator() ( MatrixType const& Af ) const
            {
                ResultType result;

                typename CoarseningType::ConnectionType connection = CoarseningType ( this->coarsening_settings_ )( result, Af );
                GridTransferType ( this->grid_transfer_settings_ )( result, connection, Af );
                CoarseMatrixConstructionType ( this->coarse_matrix_construction_settings_ )( result, Af );

                return result;
            }

          private:

            /// Settings of the coarsening
            CoarseningSettings coarsening_settings_;

            /// Settings of the interpolation
            GridTransferSettings grid_transfer_settings_;

            /// Settings of the coarse matrix construction
            CoarseMatrixConstructionSettings coarse_matrix_construction_settings_;

        };

    } // namespace AMG
} // namespace hiflow

#endif /* SRC_LINEAR_SOLVER_AMG_LEVELGENERATOR_H_ */
