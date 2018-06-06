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

#include "hiflow.h"

#ifndef MET_FLOW_INCOMP_CYL_DUAL_ASSEMBLER_H
#    define MET_FLOW_INCOMP_CYL_DUAL_ASSEMBLER_H

///
/// \brief Assembler class for incompressible Navier Stokes equations in cylindrical coordinates
///
/// \author Philipp Gerstner
///

#    include <cmath>
#    include <utility>
#    include <string>
#    include <vector>
#    include <string>
#    include <mpi.h>
#    include <sstream>
#    include <algorithm>

#    include "hiflow.h"
#    include "../../../general/goal_functional.h"
#    include "../../met_flow_incomp_assembler.h"

template<int DIM, class DataType>
class MetFlowIncompCylDualAssembler : public virtual MetFlowIncompAssembler<DIM, DataType>
{
  public:

    MetFlowIncompCylDualAssembler ( );

    ~MetFlowIncompCylDualAssembler ( )
    {
        ;
    }

    // Cellwise matrix assembly

    virtual void operator() ( const Element<double>& element, const Quadrature<double>& quadrature, LocalMatrix& lm )
    {
        MetFlowIncompAssembler<DIM, DataType>::initialize_for_element ( element, quadrature );
        MetFlowIncompCylDualAssembler<DIM, DataType>::assemble_local_matrix ( element, lm );
    }

    // cellwise vector assembly

    virtual void operator() ( const Element<double>& element, const Quadrature<double>& quadrature, LocalVector& lv )
    {
        MetFlowIncompAssembler<DIM, DataType>::initialize_for_element ( element, quadrature );
        MetFlowIncompCylDualAssembler<DIM, DataType>::assemble_local_vector ( element, lv );
    }

    virtual void operator() ( const Element<double>& element, const Quadrature<double>& quadrature, double& ls )
    {
        MetFlowIncompAssembler<DIM, DataType>::initialize_for_element ( element, quadrature );
        MetFlowIncompCylDualAssembler<DIM, DataType>::assemble_local_scalar ( element, ls );
    }

    virtual void operator() ( const Element<DataType>& element, int facet_number, const Quadrature<DataType>& quadrature, LocalVector& lv )
    {
    }

    virtual void operator() ( const Element<DataType>& element, const Quadrature<DataType>& quadrature, LocalVector& s_tau, LocalVector& s_h )
    {
    }

  protected:
    virtual void assemble_local_matrix ( const Element<double>& element, LocalMatrix& lm ) const;
    virtual void assemble_local_vector ( const Element<double>& element, LocalVector& lv ) const;
    virtual void assemble_local_scalar ( const Element<double>& element, double& ls ) const;

    virtual void assemble_local_scalar_boundary ( const Element<DataType>& element, int facet_number, LocalVector& lv ) const
    {
    }

    ///// Specific incompressible flow functions in cylindrical /////
    virtual void assemble_local_matrix_dual ( const Element<double>& element, LocalMatrix& lm ) const;
    virtual void assemble_local_vector_dual ( const Element<double>& element, LocalVector& lv ) const;

    using MetFlowAssembler<DIM, DataType>::num_dofs;
    using MetFlowAssembler<DIM, DataType>::dof_index;
    using MetFlowAssembler<DIM, DataType>::phi;
    using MetFlowAssembler<DIM, DataType>::grad_phi;
};

#endif
