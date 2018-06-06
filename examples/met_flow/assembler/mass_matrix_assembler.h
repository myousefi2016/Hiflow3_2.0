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

#ifndef MASS_MATRIX_ASSEMBLER_H
#    define MASS_MATRIX_ASSEMBLER_H

///
/// \file met_flow.h
/// \brief Assembler base classes mass matrix
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
#    include "../tmp_config/met_flow_vars.h"

///
/// \brief Abstract base class for met flow assembler implementations.
///

template<int DIM, int VARDIM, class DataType>
class MassMatrixAssembler : public AssemblyAssistant<DIM, DataType>
{
  public:
    MassMatrixAssembler ( );

    ~MassMatrixAssembler ( )
    {
    }

    virtual void operator() ( const Element<double>& element, const Quadrature<double>& quadrature, LocalMatrix& lm )
    {
        this->initialize_for_element ( element, quadrature );
        this->assemble_local_matrix ( element, lm );
    }

    virtual void initialize_for_element ( const Element<double>& element, const Quadrature<double>& quadrature )
    {
        AssemblyAssistant<DIM, DataType>::initialize_for_element ( element, quadrature );
    }
    virtual void assemble_local_matrix ( const Element<double>& element, LocalMatrix& lm ) const = 0;

  protected:

};

#endif
