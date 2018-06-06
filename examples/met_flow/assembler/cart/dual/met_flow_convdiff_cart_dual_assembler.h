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

#ifndef MET_FLOW_CONVDIFF_CART_DUAL_ASSEMBLER_H
#    define MET_FLOW_CONVDIFF_CART_DUAL_ASSEMBLER_H

///
/// \brief Assembler base class for dual scalar convection-diffusion-reaction equations in cartesian coordinates
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
#    include "../../met_flow_convdiff_assembler.h"

template<int DIM, class DataType>
class MetFlowConvDiffCartDualAssembler : public MetFlowConvDiffAssembler<DIM, DataType>
{
  public:

    MetFlowConvDiffCartDualAssembler ( );

    ~MetFlowConvDiffCartDualAssembler ( )
    {
        ;
    }

    //////////////////////////////////////////
    //// Functions from MetFlowAssembler /////
    // Cellwise matrix assembly

    virtual void operator() ( const Element<DataType>& element, const Quadrature<DataType>& quadrature, LocalMatrix& lm )
    {
        MetFlowConvDiffAssembler<DIM, DataType>::initialize_for_element ( element, quadrature );
        MetFlowConvDiffCartDualAssembler<DIM, DataType>::assemble_local_matrix ( element, lm );
    }

    // cellwise vector assembly

    virtual void operator() ( const Element<DataType>& element, const Quadrature<DataType>& quadrature, LocalVector& lv )
    {
        MetFlowConvDiffAssembler<DIM, DataType>::initialize_for_element ( element, quadrature );
        MetFlowConvDiffCartDualAssembler<DIM, DataType>::assemble_local_vector ( element, lv );
    }

    // Cellwise scalar assembly

    virtual void operator() ( const Element<DataType>& element, const Quadrature<DataType>& quadrature, DataType& ls )
    {
    }

    virtual void operator() ( const Element<DataType>& element, int facet_number, const Quadrature<DataType>& quadrature, LocalVector& lv )
    {
    }

    virtual void operator() ( const Element<DataType>& element, const Quadrature<DataType>& quadrature, LocalVector& s_tau, LocalVector& s_h )
    {
    }

    virtual void assemble_local_matrix ( const Element<DataType>& element, LocalMatrix& lm ) const;
    virtual void assemble_local_vector ( const Element<DataType>& element, LocalVector& lv ) const;

    virtual void assemble_local_scalar ( const Element<DataType>& element, DataType& ls ) const
    {
        ;
    }

    virtual void assemble_local_scalar_boundary ( const Element<DataType>& element, int facet_number, LocalVector& lv ) const
    {
        ;
    }

  protected:

    ///// Specific incompressible flow functions /////
    virtual void assemble_local_matrix_dual ( const Element<DataType>& element, LocalMatrix& lm ) const;
    virtual void assemble_local_vector_dual ( const Element<DataType>& element, LocalVector& lv ) const;

    using MetFlowAssembler<DIM, DataType>::num_dofs;
    using MetFlowAssembler<DIM, DataType>::dof_index;
    using MetFlowAssembler<DIM, DataType>::phi;
    using MetFlowAssembler<DIM, DataType>::grad_phi;
};

#endif
