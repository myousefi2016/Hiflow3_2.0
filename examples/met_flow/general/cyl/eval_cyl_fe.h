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

#ifndef EVAL_CYL_FE_H
#    define EVAL_CYL_FE_H

///
/// \file eval_cyl_fe.h
/// \brief Function for evaluating a FE function in cylindrical coordinates
///
/// \author Martin Baumann, Philipp Gerstner
///

#    include <cmath>
#    include <utility>
#    include <string>
#    include <vector>
#    include <mpi.h>

#    include "hiflow.h"
#    include "../met_flow_main.h"

using namespace hiflow;
using namespace hiflow::mesh;
using namespace hiflow::doffem;
using namespace hiflow::la;

class EvalCylFeFunction
{
  public:

    EvalCylFeFunction ( const VectorSpace<DATATYPE>& space,
                        LAD::VectorType & sol, //const PpVector<LAD>& fun,
                        int var, int rank )
    : space_ ( space ), var_ ( var ), sol_ ( sol ), rank_ ( rank )
    {
        assert ( var_ == 0 | var_ == 1 );

        //   fun_.GetDofsAndValues(all_ids,all_values);
        // std::cout << *( std::max_element( all_values.begin(), all_values.end() ) )<< "  "  << *( std::min_element( all_values.begin(), all_values.end() ) ) <<"\n";

    }

    void operator() ( const Entity& cell,
            const std::vector<double>& ref_coords,
            std::vector<double>& values ) const;

  private:
    std::vector<double> all_values;
    std::vector<int> all_ids;
    const VectorSpace<DATATYPE>& space_;
    LAD::VectorType & sol_;
    int var_;
    int rank_;
};

#endif
