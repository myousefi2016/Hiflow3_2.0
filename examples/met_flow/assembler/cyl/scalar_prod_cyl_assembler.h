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

#ifndef SCALAR_PROD_CYL_ASSEMBLERS_H
#    define SCALAR_PROD_CYL_ASSEMBLERS_H

///
/// \file met_flow.h
/// \brief Assembler base classes for L2, H1, H2 inner products
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
#    include "../scalar_prod_assembler.h"

///
/// \brief Abstract base class for met flow assembler implementations.
///

template<int DIM, int VARDIM, class DataType>
class ScalarProdCylAssembler : public ScalarProdAssembler<DIM, VARDIM, DataType>
{
  public:
    ScalarProdCylAssembler ( );

    ~ScalarProdCylAssembler ( )
    {
    }

    virtual void assemble_local_scalar ( const Element<DataType>& element, DataType& ls ) const;

  protected:

};

#endif
