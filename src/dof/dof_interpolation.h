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

#ifndef _DOF_DOF_INTERPOLATION_H_
#    define _DOF_DOF_INTERPOLATION_H_

#    include <vector>
#    include <utility>
#    include <map>
//#include <tr1/unordered_map>

#    include "common/macros.h"
#    include "dof/dof_fem_types.h"

namespace hiflow
{
    namespace doffem
    {

        ///
        /// \class  DofInterpolation dof_interpolation.h
        /// \brief  Interpolation of DoF-IDs due to hanging nodes (h- or p-refinement).
        /// \author Michael Schick<br>Martin Baumann
        ///
        /// This class is used to represent the Interpolation and Identification 
        /// between DoFs. One application area is to store (global) DoF IDs for the
        /// complete (sub-)domain another application area is to represent interpolation
        /// and identification for the interface patterns that are used during the
        /// DoF numbering procedure DegreeOfFreedom::identify_common_dofs().
        ///

        typedef std::map<DofID, std::vector<std::pair<DofID, double> > > InterpolationData;

        class DofInterpolation : public InterpolationData
        {
          public:

            typedef InterpolationData::iterator iterator;
            typedef InterpolationData::reverse_iterator reverse_iterator;
            typedef InterpolationData::const_iterator const_iterator;

            /// add interpolation definition: val(id) = sum_i val(dofs(i))*weights(i)
            bool push_back ( DofID id,
                             std::vector<DofID> const& dofs,
                             std::vector<double> const& weights );

            /// add interpolation definition
            bool push_back ( std::pair<DofID, std::vector<std::pair<DofID, double> > > interpolation );

            /// apply permutation of DoF IDs
            void apply_permutation ( const std::vector<int>& permutation );

            /// clears interpolation and identification information
            void clear_entries ( );

            /// save interpolation data to file
            void backup ( const std::string& ) const;
            void backup ( std::ostream& ) const;

            /// restore interpolation data from file
            void restore ( const std::string& );
            void restore ( std::istream& );

            std::vector<std::pair<DofID, DofID> > const& dof_identification_list ( ) const
            {
                return dof_identification_list_;
            }

            /// Fill list of dof identifications

            void insert_dof_identification ( DofID dof_slave, DofID dof_master )
            {
                dof_identification_list_.push_back ( std::make_pair ( dof_slave, dof_master ) );
            }

            bool is_constrained ( DofID node ) const;

            void add_dependency ( DofID node, DofID dep, double weight );
            void remove_dependency ( DofID node, DofID dep );
            bool has_dependency ( DofID node, DofID dep ) const;
            double get_weight ( DofID node, DofID dep ) const;
            double sum_weights ( DofID node ) const;

            void reduce ( );

            /// overloaded out stream operator
            friend std::ostream& operator<< ( std::ostream&, const DofInterpolation& );

          private:
            void reduce_from_above ( );
            void reduce_from_below ( );

            /// list contains corresponding DoFs that should be identified
            std::vector<std::pair<DofID, DofID> > dof_identification_list_;

        };

    }
} // namespace hiflow
#endif
