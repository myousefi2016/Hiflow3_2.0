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

#include "mass_matrix_assembler.h"

template<int DIM, int VARDIM, class DataType>
MassMatrixAssembler<DIM, VARDIM, DataType>::MassMatrixAssembler ( )
{
}

/// Template instanciation.
template class MassMatrixAssembler<2, 1, double>;
template class MassMatrixAssembler<3, 1, double>;
template class MassMatrixAssembler<2, 2, double>;
template class MassMatrixAssembler<3, 2, double>;
template class MassMatrixAssembler<2, 3, double>;
template class MassMatrixAssembler<3, 3, double>;
template class MassMatrixAssembler<2, 4, double>;
template class MassMatrixAssembler<3, 4, double>;
template class MassMatrixAssembler<2, 5, double>;
template class MassMatrixAssembler<3, 5, double>;
template class MassMatrixAssembler<2, 6, double>;
template class MassMatrixAssembler<3, 6, double>;
template class MassMatrixAssembler<2, 7, double>;
template class MassMatrixAssembler<3, 7, double>;
