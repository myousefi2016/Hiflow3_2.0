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

#include "met_flow_incomp_est_assembler.h"

template<int DIM, class DataType>
MetFlowIncompEstimatorAssembler<DIM, DataType>::MetFlowIncompEstimatorAssembler ( )
: MetFlowEstimatorAssembler<DIM, DataType>( ),
MetFlowIncompAssembler<DIM, DataType>( )
{
    this->mode_ = ESTIMATOR;
    this->num_eq_ = DIM + 1;
#ifdef AUGMENT_PRESS
    this->num_eq_++;
#endif

    MetFlowEstimatorAssembler<DIM, DataType>::num_var_ = DIM + 1;
#ifdef AUGMENT_PRESS
    MetFlowEstimatorAssembler<DIM, DataType>::num_var_++;
#endif
    const int num_var = MetFlowEstimatorAssembler<DIM, DataType>::num_var_;

    MetFlowEstimatorAssembler<DIM, DataType>::allocate_function_evaluators ( num_var );
    MetFlowEstimatorAssembler<DIM, DataType>::allocate_function_values ( num_var );
}

template<int DIM, class DataType>
void MetFlowIncompEstimatorAssembler<DIM, DataType>::clear ( )
{
    MetFlowEstimatorAssembler<DIM, DataType>::clear ( );
    MetFlowIncompAssembler<DIM, DataType>::clear ( );

    this->mode_ = ESTIMATOR;
    this->num_eq_ = DIM + 1;
#ifdef AUGMENT_PRESS
    this->num_eq_++;
#endif

    MetFlowEstimatorAssembler<DIM, DataType>::num_var_ = DIM + 1;
#ifdef AUGMENT_PRESS
    MetFlowEstimatorAssembler<DIM, DataType>::num_var_++;
#endif
    const int num_var = MetFlowEstimatorAssembler<DIM, DataType>::num_var_;

    MetFlowEstimatorAssembler<DIM, DataType>::allocate_function_evaluators ( num_var );
    MetFlowEstimatorAssembler<DIM, DataType>::allocate_function_values ( num_var );
}

template class MetFlowIncompEstimatorAssembler<2, double>;
template class MetFlowIncompEstimatorAssembler<3, double>;

