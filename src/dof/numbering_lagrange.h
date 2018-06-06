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

#ifndef _DOF_NUMBERING_LAGRANGE_H_
#    define _DOF_NUMBERING_LAGRANGE_H_

#    include "numbering_strategy.h"

namespace hiflow
{
    namespace doffem
    {

        /// \author Michael Schick<br>Martin Baumann

        template<class DataType>
        class NumberingLagrange : public NumberingStrategy<DataType>
        {
          public:

            /// Dense matrix type for interpolation weights.
            typedef std::vector< std::vector<DataType> > InterpolationWeights;

            /// Sparse matrix type for describing interpolation constraints on
            /// subset of DoF:s lying on an interface.
            typedef std::map<int, std::vector<DataType> > InterfaceMatrix;

            /// Constructor

            NumberingLagrange ( ) : NumberingStrategy<DataType>( )
            {
            }

            /// Implementation of numbering procedure for the given mesh
            /// @param[in] order Ordering strategy for DoFs.
            void number ( DOF_ORDERING order = HIFLOW_CLASSIC );

          private:

            /// Discontinuous initial numbering of all occuring dofs within the given mesh
            void initial_numbering ( );
            /// Dofs are idendified due to constraints given by the Lagrange FE Type
            void identify_common_dofs ( );

            /// Status information about interface patterns
            void print_interface_patterns ( ) const;

            /// Mapping needed for identification of dofs
            std::map<FEInterfacePattern<DataType>, DofInterpolationPattern>
            interface_patterns_interpolation_;

            /// Compute the interpolation weights, which are needed for identifications
            void compute_interpolation ( const FEInterfacePattern<DataType>& pattern,
                                         const mesh::Interface& interface,
                                         DofInterpolationPattern& interpolation );

            /// Calculate weights of DoFs for two neighbouring cells, i.e. B-DoFs are 
            /// represented by weighted sum of A-DoFs
            void compute_weights ( mesh::Entity const& cellA,
                                   doffem::FEType<DataType> const& ansatzA,
                                   mesh::Entity const& cellB,
                                   doffem::FEType<DataType> const& ansatzB,
                                   InterpolationWeights& weights ) const;

        };

    }
}

#endif
