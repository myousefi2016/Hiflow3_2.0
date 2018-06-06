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

#include "lmatrix.h"
#include "lmatrix_csr.h"

// class CSR_lMatrix

template <typename ValueType>
CSR_lMatrix<ValueType>::CSR_lMatrix ( )
{
    this->matrix_format_name_ = "CSR matrix";
    this->matrix_format_id_ = hiflow::la::CSR;
}

template <typename ValueType>
CSR_lMatrix<ValueType>::~CSR_lMatrix ( )
{
}

template class CSR_lMatrix<float>;
template class CSR_lMatrix<double>;
