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

#include "met_flow_convdiff_est_assembler.h"

template<int DIM, class DataType>
MetFlowConvDiffEstimatorAssembler<DIM, DataType>::MetFlowConvDiffEstimatorAssembler ( )
: MetFlowEstimatorAssembler<DIM, DataType>( ),
MetFlowConvDiffAssembler<DIM, DataType>( )
{
    this->num_eq_ = 1;
    MetFlowEstimatorAssembler<DIM, DataType>::num_var_ = 1;
    const int num_var = MetFlowEstimatorAssembler<DIM, DataType>::num_var_;

    MetFlowEstimatorAssembler<DIM, DataType>::allocate_function_evaluators ( num_var );
    MetFlowEstimatorAssembler<DIM, DataType>::allocate_function_values ( num_var );
}

template<int DIM, class DataType>
void MetFlowConvDiffEstimatorAssembler<DIM, DataType>::clear ( )
{
    MetFlowEstimatorAssembler<DIM, DataType>::clear ( );
    MetFlowConvDiffAssembler<DIM, DataType>::clear ( );

    this->num_eq_ = 1;
    MetFlowEstimatorAssembler<DIM, DataType>::num_var_ = 1;
    const int num_var = MetFlowEstimatorAssembler<DIM, DataType>::num_var_;

    MetFlowEstimatorAssembler<DIM, DataType>::allocate_function_evaluators ( num_var );
    MetFlowEstimatorAssembler<DIM, DataType>::allocate_function_values ( num_var );
}

template class MetFlowConvDiffEstimatorAssembler<2, double>;
template class MetFlowConvDiffEstimatorAssembler<3, double>;

