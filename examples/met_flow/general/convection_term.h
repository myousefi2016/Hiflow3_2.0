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

#ifndef CONVECTION_TERM_H_
#    define CONVECTION_TERM_H_

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
class ConvectionTerm
{
  public:

    ConvectionTerm ( )
    {
    }

    virtual void evaluate ( const DataType t, const Vec<DIM, DataType>& x, const int d, DataType& val ) const = 0;

    virtual void evaluate_grad ( const DataType t, const Vec<DIM, DataType>& x, const int d, Vec<DIM, DataType>& val ) const = 0;
};

template class ConvectionTerm< 2, double>;
template class ConvectionTerm< 3, double>;

template<int DIM, class DataType>
class ConstantConvectionTerm : public virtual ConvectionTerm<DIM, DataType>
{
  public:

    ConstantConvectionTerm ( Vec<DIM, DataType> const_conv )
    : conv_ ( const_conv )
    {
    }

    virtual void evaluate ( const DataType t, const Vec<DIM, DataType>& x, const int d, DataType& val ) const
    {
        val = conv_[d];
    }

    virtual void evaluate_grad ( const DataType t, const Vec<DIM, DataType>& x, const int d, Vec<DIM, DataType>& val ) const
    {
        val = Vec<DIM, DataType>( );
    }

    const Vec<DIM, DataType> conv_;
};

template class ConstantConvectionTerm< 2, double>;
template class ConstantConvectionTerm< 3, double>;

template<int DIM, class DataType>
class RotationConvectionTerm : public virtual ConvectionTerm<DIM, DataType>
{
  public:

    RotationConvectionTerm ( Vec<3, DataType> normal, Vec<3, DataType> center )
    : normal_ ( normal ),
    center_ ( center )
    {
    }

    virtual void evaluate ( const DataType t, const Vec<DIM, DataType>& x, const int d, DataType& val ) const
    {
        std::vector<DataType> r ( 3, 0. );
        r[0] = x[0] - center_[0];
        r[1] = x[1] - center_[1];
        if ( DIM > 2 )
        {
            r[2] = x[2] - center_[2];
        }

        switch ( d )
        {
            case 0:
                val = normal_[1] * r[2] - normal_[2] * r[1];
                break;
            case 1:
                val = normal_[2] * r[0] - normal_[0] * r[2];
                break;
            case 2:
                val = normal_[0] * r[1] - normal_[1] * r[0];
                break;
        }
    }

    virtual void evaluate_grad ( const DataType t, const Vec<DIM, DataType>& x, const int d, Vec<DIM, DataType>& val ) const
    {
        switch ( d )
        {
            case 0:
                val[0] = 0.;
                val[1] = -normal_[2];
                if ( DIM > 2 )
                {
                    val[2] = normal_[1];
                }
                break;
            case 1:
                val[0] = normal_[2];
                val[1] = 0.;
                if ( DIM > 2 )
                {
                    val[2] = -normal_[0];
                }
                break;
            case 2:
                val[0] = -normal_[1];
                val[1] = normal_[0];
                if ( DIM > 2 )
                {
                    val[2] = 0.;
                }
                break;
        }
    }

    const Vec<3, DataType> normal_;
    const Vec<3, DataType> center_;

};

template class RotationConvectionTerm< 2, double>;
template class RotationConvectionTerm< 3, double>;

#endif
