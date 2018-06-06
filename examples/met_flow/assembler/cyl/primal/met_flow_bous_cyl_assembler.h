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

#ifndef MET_FLOW_BOUS_CYL_ASSEMBLER_H
#    define MET_FLOW_BOUS_CYL_ASSEMBLER_H

///
/// \brief Assembler class for Boussinesq equations in cylindrical coordinates
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
#    include "../../met_flow_bous_assembler.h"
#    include "met_flow_incomp_cyl_assembler.h"

template<int DIM, class DataType>
class MetFlowBousCylAssembler : public virtual MetFlowBousAssembler<DIM, DataType>, public virtual MetFlowIncompCylAssembler<DIM, DataType>
{
  public:

    MetFlowBousCylAssembler ( );

    ~MetFlowBousCylAssembler ( )
    {
        ;
    }

    // Cellwise matrix assembly

    virtual void operator() ( const Element<DataType>& element, const Quadrature<DataType>& quadrature, LocalMatrix& lm )
    {
        MetFlowBousAssembler<DIM, DataType>::initialize_for_element ( element, quadrature );
        MetFlowBousCylAssembler<DIM, DataType>::assemble_local_matrix ( element, lm );
    }

    // cellwise vector assembly

    virtual void operator() ( const Element<DataType>& element, const Quadrature<DataType>& quadrature, LocalVector& lv )
    {
        MetFlowBousAssembler<DIM, DataType>::initialize_for_element ( element, quadrature );
        MetFlowBousCylAssembler<DIM, DataType>::assemble_local_vector ( element, lv );
    }

    // Cellwise scalar assembly

    virtual void operator() ( const Element<DataType>& element, const Quadrature<DataType>& quadrature, DataType& ls )
    {
        MetFlowBousAssembler<DIM, DataType>::initialize_for_element ( element, quadrature );
        MetFlowBousCylAssembler<DIM, DataType>::assemble_local_scalar ( element, ls );
    }

    virtual void operator() ( const Element<DataType>& element, int facet_number, const Quadrature<DataType>& quadrature, LocalVector& lv )
    {
        MetFlowBousAssembler<DIM, DataType>::initialize_for_facet ( element, quadrature, facet_number );
        MetFlowBousCylAssembler<DIM, DataType>::assemble_local_scalar_boundary ( element, facet_number, lv );
    }

    virtual void operator() ( const Element<DataType>& element, const Quadrature<DataType>& quadrature, LocalVector& s_tau, LocalVector& s_h )
    {
    }

    virtual void assemble_local_matrix ( const Element<DataType>& element, LocalMatrix& lm ) const;
    virtual void assemble_local_vector ( const Element<DataType>& element, LocalVector& lv ) const;
    virtual void assemble_local_scalar ( const Element<DataType>& element, DataType& ls ) const;
    virtual void assemble_local_scalar_boundary ( const Element<DataType>& element, int facet_number, LocalVector& lv ) const;

  protected:

    ////////////////////////////////////////////////
    //// Functions from MetFlowIncompAssembler /////

    virtual void assemble_local_matrix_primal ( const Element<DataType>& element, LocalMatrix& lm ) const;

    virtual void assemble_local_vector_primal ( const Element<DataType>& element, LocalVector& lv ) const;
    virtual void assemble_local_vector_quant ( const Element<DataType>& element, LocalVector& lv ) const;
    virtual void assemble_local_vector_goal ( const Element<DataType>& element, LocalVector& lv ) const;

    // TODO
    virtual void assemble_local_scalar_goal_int ( const Element<DataType>& element, DataType& ls ) const;
    virtual void assemble_local_scalar_goal_fin ( const Element<DataType>& element, DataType& ls ) const;

    ////////////////////////////////////////////////
    //// Specific functions for Boussinesq /////////
    virtual void assemble_local_vector_grad_temp ( const Element<DataType>& element, LocalVector& lv ) const;
    virtual void assemble_local_vector_buoyancy ( const Element<DataType>& element, LocalVector& lv ) const;

    // TODO
    virtual void assemble_local_scalar_temp_mean ( const Element<DataType>& element, DataType& ls ) const;
    virtual void assemble_local_scalar_rad_heat_surface ( const Element<DataType>& element, int facet_number, LocalVector& lv ) const;

    using MetFlowAssembler<DIM, DataType>::num_dofs;
    using MetFlowAssembler<DIM, DataType>::dof_index;
    using MetFlowAssembler<DIM, DataType>::phi;
    using MetFlowAssembler<DIM, DataType>::grad_phi;
};

#endif
