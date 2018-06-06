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

#ifndef MET_FLOW_CONVDIFF_CART_ESTIMATOR_ASSEMBLER_H
#    define MET_FLOW_CONVDIFF_CART_ESTIMATOR_ASSEMBLER_H

///
/// \brief Assembler class for DWR estimator for convection diffusion equation in cartesian coordinates
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
#    include "../../estimator/met_flow_convdiff_est_assembler.h"

template<int DIM, class DataType>
class MetFlowConvDiffCartEstimatorAssembler : public virtual MetFlowConvDiffEstimatorAssembler<DIM, DataType>
{
    typedef typename DGGlobalAssembler<DataType>::InterfaceSide InterfaceSide;

  public:

    MetFlowConvDiffCartEstimatorAssembler ( );

    ~MetFlowConvDiffCartEstimatorAssembler ( )
    {
        ;
    }

  protected:

    virtual void assemble_primal_cell_residual ( const Element<DataType>& element, int eq, int est_mode, DataType& r, DataType& w ) const;
    virtual void assemble_dual_cell_residual ( const Element<DataType>& element, int eq, int est_mode, DataType& r, DataType& w ) const;
    virtual void assemble_dual_cell_time_jump ( const Element<DataType>& element, int eq, int est_mode, DataType& r, DataType& w ) const;

    virtual void assemble_primal_interface_jump ( const Element<DataType>& left_elem,
                                                  const Element<DataType>& right_elem,
                                                  const Quadrature<DataType>& left_quad,
                                                  const Quadrature<DataType>& right_quad,
                                                  int left_facet_number,
                                                  int right_facet_number,
                                                  InterfaceSide left_if_side,
                                                  InterfaceSide right_if_side,
                                                  int eq,
                                                  int est_mode,
                                                  DataType& r,
                                                  DataType& w );

    virtual void assemble_dual_interface_jump ( const Element<DataType>& left_elem,
                                                const Element<DataType>& right_elem,
                                                const Quadrature<DataType>& left_quad,
                                                const Quadrature<DataType>& right_quad,
                                                int left_facet_number,
                                                int right_facet_number,
                                                InterfaceSide left_if_side,
                                                InterfaceSide right_if_side,
                                                int eq,
                                                int est_mode,
                                                DataType& r,
                                                DataType& w );

    virtual void assemble_primal_interface_boundary ( const Element<DataType>& left_elem,
                                                      const Element<DataType>& right_elem,
                                                      const Quadrature<DataType>& left_quad,
                                                      const Quadrature<DataType>& right_quad,
                                                      int left_facet_number,
                                                      int right_facet_number,
                                                      InterfaceSide left_if_side,
                                                      InterfaceSide right_if_side,
                                                      int eq,
                                                      int est_mode,
                                                      DataType& r,
                                                      DataType& w );

    virtual void assemble_dual_interface_boundary ( const Element<DataType>& left_elem,
                                                    const Element<DataType>& right_elem,
                                                    const Quadrature<DataType>& left_quad,
                                                    const Quadrature<DataType>& right_quad,
                                                    int left_facet_number,
                                                    int right_facet_number,
                                                    InterfaceSide left_if_side,
                                                    InterfaceSide right_if_side,
                                                    int eq,
                                                    int est_mode,
                                                    DataType& r,
                                                    DataType& w );

};

#endif
