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

#ifndef _PERIODIC_BOUNDARY_CONDITIONS_CARTESIAN_H_
#    define _PERIODIC_BOUNDARY_CONDITIONS_CARTESIAN_H_

#    include <string>
#    include <vector>
#    include <map>
#    include <mpi.h>

#    include "space/periodic_boundary_conditions.h"
#    include "dof/degree_of_freedom.h"

namespace hiflow
{

    ///
    /// @class PeriodicBoundaryConditionsCartesian periodic_boundary_conditions_cartesian.h
    /// @brief Handling of periodic boundary conditions in case of boundaries orthogonal to cartesian axes.
    /// @author Martin Baumann
    ///
    /// Only the special case is considered where periodic boundaries are orthogonal
    /// to x, y or z axis. Therefore the matching DoFs can be identified by looking
    /// at the DoFs coordinates.
    ///

    template<class DataType>
    class PeriodicBoundaryConditionsCartesian : public PeriodicBoundaryConditions<DataType>
    {
      protected:

        std::map<doffem::DofID, DataType> x_coord_;
        std::map<doffem::DofID, DataType> y_coord_;
        std::map<doffem::DofID, DataType> z_coord_;

        void compute_conditions ( VectorSpace<DataType>& space );

      public:

        PeriodicBoundaryConditionsCartesian ( )
        {
        };

        ~PeriodicBoundaryConditionsCartesian ( )
        {
        };

    };

} // namespace hiflow

#endif
