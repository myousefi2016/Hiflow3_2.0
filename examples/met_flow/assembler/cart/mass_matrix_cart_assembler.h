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

#ifndef MASS_MATRIX_CART_ASSEMBLER_H
#    define MASS_MATRIX_CART_ASSEMBLER_H

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
#    include "../mass_matrix_assembler.h"

template<int DIM, int VARDIM, class DataType>
class MassMatrixCartAssembler : public MassMatrixAssembler<DIM, VARDIM, DataType>
{
  public:

    MassMatrixCartAssembler ( )
    {
    }

    ~MassMatrixCartAssembler ( )
    {
    }

    virtual void assemble_local_matrix ( const Element<DataType>& element, LocalMatrix& lm ) const;
  protected:

};

#endif
