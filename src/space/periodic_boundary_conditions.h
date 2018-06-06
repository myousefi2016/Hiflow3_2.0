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

#ifndef _PERIODIC_BOUNDARY_CONDITIONS_H_
#    define _PERIODIC_BOUNDARY_CONDITIONS_H_

#    include <string>
#    include <vector>
#    include <map>
#    include <mpi.h>

#    include "dof/degree_of_freedom.h"
#    include "dof/dof_fem_types.h"
#    include "space/vector_space.h"

namespace hiflow
{

    ///
    /// @class PeriodicBoundaryConditions periodic_boundary_conditions.h
    /// @brief Handling of periodic boundary conditions.
    /// @author Martin Baumann
    ///

    /**
    Periodic boundary conditions are applied by identification of degrees of freedom,
    hence no interpolation or iterative procedure is needed.
    <H2>Example:</H2>
    <img src="../images/periodic_boundary.png" alt="DoF identification in case of periodic boundaries">
    This class serves as abstract class, as the member function compute_conditions()
    is pure abstract. This function calculates the corresponding degrees of freedom.
    By the material number of boundaries, the periodic boundary condition can be set.
    A doubly periodic domain can be set by adding two boundary tuples, f.e..
    The colours of the DoF points on the boundaries indicate unique DoFs, i.e.
    the green DoFs on the left-hand side are unified with the ones on the
    right-hand side. Due to the doubly periodicity condition, the four blue DoFs
    are unified to one single DoF.
    \code
      VectorSpace space;
      // ...
      PeriodicBoundaryConditionsCartesian periodic_boundaries;
      periodic_boundaries.add_boundary_tuple(10, 12);
      periodic_boundaries.add_boundary_tuple(11, 13);
      periodic_boundaries.apply_boundary_conditions(space);
      std::cout << periodic_boundaries << std::endl;
    \endcode
    \todo Should there be a base class for general boundary condition
          that can iterate over boundaries of the real domain?
     **/

    template<class DataType>
    class PeriodicBoundaryConditions
    {
      public:

        typedef mesh::MaterialNumber MaterialNumber;

      protected:

        /// maps a 'real' dof to vector of corresponding dofs, which will be eliminated
        std::map<doffem::DofID, std::vector<doffem::DofID> > corresponding_dof_;

        /// mapping (mat,var) -> dof_id
        std::map<std::pair<MaterialNumber, int>, std::vector<doffem::DofID> > dof_id_map_;

        /// find corresponding DoFs if there are more than one
        void handle_multiple_correspondences ( );

        /// description of boundary correspondence, f.e. mat 12 <-> mat 13
        std::map<MaterialNumber, MaterialNumber> boundary_descriptor_;

        /// returns true if entry was added; returns false if entry already existed.
        bool insert_dof_id_entry ( MaterialNumber mat, int var, doffem::DofID dof_id );

        /// fills the list of corresponding DoFs (corresponding_dof_)
        virtual void compute_conditions ( VectorSpace<DataType>& ) = 0;

        /// permutes DoFs such that number of DoFs is reduced
        void change_dofs ( VectorSpace<DataType>& space );

      public:

        PeriodicBoundaryConditions ( )
        {
        };

        virtual ~PeriodicBoundaryConditions ( )
        {
        };

        /// add boundary tuple to boundary_descriptor_
        void add_boundary_tuple ( MaterialNumber, MaterialNumber );

        /// calculates corresponding DoFs and performs DoF identification
        void apply_boundary_conditions ( VectorSpace<DataType>& );

        /// checks if a boundary has periodic boundary condition
        bool is_periodic_boundary ( MaterialNumber ) const;

    };

} // namespace hiflow

#endif
