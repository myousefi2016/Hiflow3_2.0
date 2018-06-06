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

#ifndef MET_FLOW_INCOMP_CYL_ASSEMBLER_H
#    define MET_FLOW_INCOMP_CYL_ASSEMBLER_H

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
#    include "../../met_flow_incomp_assembler.h"

template<int DIM, class DataType>
class MetFlowIncompCylAssembler : public virtual MetFlowIncompAssembler<DIM, DataType>
{
  public:

    MetFlowIncompCylAssembler ( );

    ~MetFlowIncompCylAssembler ( )
    {
        ;
    }

    // Cellwise matrix assembly

    virtual void operator() ( const Element<DataType>& element, const Quadrature<DataType>& quadrature, LocalMatrix& lm )
    {
        MetFlowIncompAssembler<DIM, DataType>::initialize_for_element ( element, quadrature );
        MetFlowIncompCylAssembler<DIM, DataType>::assemble_local_matrix ( element, lm );
    }

    // cellwise vector assembly

    virtual void operator() ( const Element<DataType>& element, const Quadrature<DataType>& quadrature, LocalVector& lv )
    {
        MetFlowIncompAssembler<DIM, DataType>::initialize_for_element ( element, quadrature );
        MetFlowIncompCylAssembler<DIM, DataType>::assemble_local_vector ( element, lv );
    }

    virtual void operator() ( const Element<DataType>& element, const Quadrature<DataType>& quadrature, DataType& ls )
    {
        MetFlowIncompAssembler<DIM, DataType>::initialize_for_element ( element, quadrature );
        MetFlowIncompCylAssembler<DIM, DataType>::assemble_local_scalar ( element, ls );
    }

    virtual void operator() ( const Element<DataType>& element, int facet_number, const Quadrature<DataType>& quadrature, LocalVector& lv )
    {
        MetFlowIncompAssembler<DIM, DataType>::initialize_for_facet ( element, quadrature, facet_number );
        MetFlowIncompCylAssembler<DIM, DataType>::assemble_local_scalar_boundary ( element, facet_number, lv );
    }

    virtual void operator() ( const Element<DataType>& element, const Quadrature<DataType>& quadrature, LocalVector& s_tau, LocalVector& s_h )
    {
    }

    virtual void assemble_local_matrix ( const Element<DataType>& element, LocalMatrix& lm ) const;
    virtual void assemble_local_vector ( const Element<DataType>& element, LocalVector& lv ) const;
    virtual void assemble_local_scalar ( const Element<DataType>& element, DataType& ls ) const;
    virtual void assemble_local_scalar_boundary ( const Element<DataType>& element, int facet_number, LocalVector& lv ) const;

  protected:
    ///// Specific incompressible flow functions in cylindrical /////
    virtual void assemble_local_matrix_primal ( const Element<DataType>& element, LocalMatrix& lm ) const;

    virtual void assemble_local_vector_primal ( const Element<DataType>& element, LocalVector& lv ) const;
    virtual void assemble_local_vector_vort ( const Element<DataType>& element, LocalVector& lv ) const;
    virtual void assemble_local_vector_quant ( const Element<DataType>& element, LocalVector& lv ) const;
    virtual void assemble_local_vector_goal ( const Element<DataType>& element, LocalVector& lv ) const;

    virtual void assemble_local_scalar_goal_int ( const Element<DataType>& element, DataType& ls ) const;
    virtual void assemble_local_scalar_goal_fin ( const Element<DataType>& element, DataType& ls ) const;
    virtual void assemble_local_scalar_div_max ( const Element<DataType>& element, DataType& ls ) const;
    virtual void assemble_local_scalar_div_mean ( const Element<DataType>& element, DataType& ls ) const;

    using MetFlowAssembler<DIM, DataType>::num_dofs;
    using MetFlowAssembler<DIM, DataType>::dof_index;
    using MetFlowAssembler<DIM, DataType>::phi;
    using MetFlowAssembler<DIM, DataType>::grad_phi;
};

#endif
