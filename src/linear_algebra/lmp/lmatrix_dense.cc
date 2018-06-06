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

/// @author Niels Wegh

#include "lmatrix.h"
#include "lmatrix_dense.h"

// class DENSE_lMatrix

template <typename ValueType>
DENSE_lMatrix<ValueType>::DENSE_lMatrix ( )
{
    this->matrix_format_name_ = "DENSE";
    this->matrix_format_id_ = hiflow::la::DENSE;
}

template <typename ValueType>
DENSE_lMatrix<ValueType>::~DENSE_lMatrix ( )
{
}

template class DENSE_lMatrix<float>;
template class DENSE_lMatrix<double>;
