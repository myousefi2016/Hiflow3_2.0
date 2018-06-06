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

#ifndef MET_FLOW_CONVDIFF_EST_ASSEMBLER_H
#    define MET_FLOW_CONVDIFF_EST_ASSEMBLER_H

///
/// \brief Assembler base class for scalar convection-diffusion-reaction equations
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
#    include "../met_flow_convdiff_assembler.h"
#    include "../met_flow_est_assembler.h"

template<int DIM, class DataType>
class MetFlowConvDiffEstimatorAssembler : public MetFlowEstimatorAssembler<DIM, DataType>, public MetFlowConvDiffAssembler<DIM, DataType>
{
    typedef typename DGGlobalAssembler<DataType>::InterfaceSide InterfaceSide;

  public:

    MetFlowConvDiffEstimatorAssembler ( );

    ~MetFlowConvDiffEstimatorAssembler ( )
    {
        ;
    }

    virtual void operator() ( const Element<DataType>& element, const Quadrature<DataType>& quadrature, LocalVector& vals )
    {
        MetFlowEstimatorAssembler<DIM, DataType>::assemble_local_est_cell ( element, quadrature, vals );
    }

    virtual void operator() ( const Element<DataType>& left_elem,
            const Element<DataType>& right_elem,
            const Quadrature<DataType>& left_quad,
            const Quadrature<DataType>& right_quad,
            int left_facet_number,
            int right_facet_number,
            InterfaceSide left_if_side,
            InterfaceSide right_if_side,
            LocalVector& vals )
    {
        MetFlowEstimatorAssembler<DIM, DataType>::assemble_local_est_interface ( left_elem, right_elem, left_quad, right_quad, left_facet_number,
                                                                                 right_facet_number, left_if_side, right_if_side, vals );
    }

    virtual void clear ( );

  protected:

};

#endif
