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

#ifndef SCALAR_PROD_ASSEMBLER_H
#    define SCALAR_PROD_ASSEMBLER_H

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
#    include "../tmp_config/met_flow_vars.h"

///
/// \brief Abstract base class for met flow assembler implementations.
///

template<int DIM, int VARDIM, class DataType>
class ScalarProdAssembler : public AssemblyAssistant<DIM, DataType>
{
  public:
    ScalarProdAssembler ( );

    ~ScalarProdAssembler ( )
    {
    }

    void set_left_vector ( const LAD::VectorType& vec )
    {
        this->vector_left_ = &vec;
    }

    void set_right_vector ( const LAD::VectorType& vec )
    {
        this->vector_right_ = &vec;
    }

    void set_mode ( std::vector<bool>& L2, std::vector<bool>& H1, std::vector<bool>& H2 )
    {
        assert ( L2.size ( ) == VARDIM );
        assert ( H1.size ( ) == VARDIM );
        assert ( H2.size ( ) == VARDIM );
        this->L2_ = L2;
        this->H1_ = H1;
        this->H2_ = H2;
    }

    void mark_vector_field ( std::vector<int> indices )
    {
        assert ( indices.size ( ) == DIM );
        this->vector_field_ind_.push_back ( indices );
        for ( int l = 0; l < DIM; ++l )
        {
            this->belongs_to_field_[indices[l]] = true;
        }
    }

    virtual void operator() ( const Element<DataType>& element, const Quadrature<DataType>& quadrature, LocalMatrix& lm )
    {
    };

    virtual void operator() ( const Element<DataType>& element, const Quadrature<DataType>& quadrature, LocalVector& lv )
    {
    };

    virtual void operator() ( const Element<DataType>& element, const Quadrature<DataType>& quadrature, DataType& ls )
    {
        ls = 0.;
        this->initialize_for_element ( element, quadrature );
        this->assemble_local_scalar ( element, ls );
    };

    virtual void operator() ( const Element<DataType>& element, int facet_number, const Quadrature<DataType>& quadrature, DataType& ls )
    {
        ;
    }

    virtual void operator() ( const Element<DataType>& element, const Quadrature<DataType>& quadrature, LocalVector& s_tau, LocalVector& s_h )
    {
        ;
    }

    virtual void initialize_for_element ( const Element<DataType>& element, const Quadrature<DataType>& quadrature );

    virtual void initialize_for_facet ( const Element<DataType>& element, const Quadrature<DataType>& quadrature, int facet_number )
    {
        ;
    }

    virtual void initialize_for_element_estimator ( const Element<DataType>& element, const Quadrature<DataType>& quadrature )
    {
        ;
    }

    virtual void initialize_for_facet_estimator ( const Element<DataType>& element, const Quadrature<DataType>& quadrature, int facet_number )
    {
        ;
    }

    virtual void assemble_local_matrix ( const Element<DataType>& element, LocalMatrix& lm ) const
    {
        ;
    }

    virtual void assemble_local_vector ( const Element<DataType>& element, LocalVector& lv ) const
    {
        ;
    }
    virtual void assemble_local_scalar ( const Element<DataType>& element, DataType& ls ) const = 0;

    virtual void assemble_local_scalar_boundary ( const Element<DataType>& element, int facet_number, DataType& ls ) const
    {
        ;
    }

  protected:
    std::vector<bool> L2_;
    std::vector<bool> H1_;
    std::vector<bool> H2_;

    LAD::VectorType const* vector_left_;
    LAD::VectorType const* vector_right_;

    FunctionValues< DataType> left_[VARDIM];
    FunctionValues< Vec<DIM, DataType> > grad_left_[VARDIM];
    FunctionValues< Mat<DIM, DIM, DataType> > hess_left_[VARDIM];

    FunctionValues< DataType> right_[VARDIM];
    FunctionValues< Vec<DIM, DataType> > grad_right_[VARDIM];
    FunctionValues< Mat<DIM, DIM, DataType> > hess_right_[VARDIM];

    std::vector< std::vector<int> > vector_field_ind_;
    std::vector<bool> belongs_to_field_;

    LAD::VectorType const& vector_left ( ) const
    {
        return *this->vector_left_;
    }

    LAD::VectorType const& vector_right ( ) const
    {
        return *this->vector_right_;
    }

};

#endif
