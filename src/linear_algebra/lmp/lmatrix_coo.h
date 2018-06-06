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

/// @author Dimitar Lukarski

#ifndef __LMATRIX_COO_H
#    define __LMATRIX_COO_H

#    include "lmatrix.h"

/// @brief The base COO class
/// @author Dimitar Lukarski

template <typename ValueType>
class COO_lMatrix : public hiflow::la::lMatrix<ValueType>
{
  public:

    COO_lMatrix ( );
    virtual ~COO_lMatrix ( );

};

#endif
