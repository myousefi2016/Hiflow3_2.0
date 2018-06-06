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

#ifndef SOURCE_TERM_H_
#    define SOURCE_TERM_H_

#    include <iostream>
#    include <string>
#    include "assembly/assembly_assistant.h"
#    include "assembly/assembly.h"
#    include "assembly/function_values.h"
#    include "common/vector_algebra.h"
#    include "assembly/function_values.h"
#    include "linear_algebra/la_descriptor.h"
#    include "../tmp_config/met_flow_vars.h"

/// Structure for the construction of dirichlet boundaries

template<int DIM, class DataType>
class SourceTerm
{
  public:

    SourceTerm ( )
    {
    }

    virtual void evaluate ( int var, const DataType t, const Vec<DIM, DataType>& x, DataType& val ) const = 0;
};

template class SourceTerm< 2, double>;
template class SourceTerm< 3, double>;

template<int DIM, class DataType>
class LocalSourceTerm : public virtual SourceTerm<DIM, DataType>
{
  public:

    LocalSourceTerm ( int var, Vec<DIM, DataType> pos, DataType width, DataType ampl, DataType start, DataType stop )
    : var_ ( var ),
    pos_ ( pos ),
    width_ ( width ),
    ampl_ ( ampl ),
    start_time_ ( start ),
    stop_time_ ( stop )
    {
    }

    virtual void evaluate ( int var, const DataType t, const Vec<DIM, DataType>& x, DataType& val ) const
    {
        if ( var != var_ )
        {
            return;
        }

        /*
        if (norm(x-pos_) > width_)
        {
            val = 0.;
            return;
        }
         */
        if ( t < start_time_ )
        {
            val = 0.;
            return;
        }
        if ( t > stop_time_ )
        {
            val = 0.;
            return;
        }
        if ( std::abs ( x[0] - pos_[0] ) > 0.5 * width_ )
        {
            val = 0.;
            return;
        }
        if ( std::abs ( x[1] - pos_[1] ) > 0.5 * width_ )
        {
            val = 0.;
            return;
        }
        if ( DIM > 2 )
        {
            if ( std::abs ( x[2] - pos_[2] ) > 0.5 * width_ )
            {
                val = 0.;
                return;
            }
        }
        //val = ampl_ * (1. - norm(x-pos_) / width_);
        val = ampl_;
    }

    const DataType width_;
    const DataType ampl_;
    const int var_;
    const Vec<DIM, DataType> pos_;
    const DataType start_time_;
    const DataType stop_time_;
};

template class LocalSourceTerm< 2, double>;
template class LocalSourceTerm< 3, double>;

#endif
