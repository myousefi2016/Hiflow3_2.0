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

#ifndef METFLOW_GENERAL_GOAL_FUNCTIONAL_H
#    define METFLOW_GENERAL_GOAL_FUNCTIONAL_H

#    include <iostream>
#    include <string>
#    include "assembly/assembly_assistant.h"
#    include "assembly/assembly.h"
#    include "assembly/function_values.h"
#    include "common/vector_algebra.h"
#    include "assembly/function_values.h"
#    include "linear_algebra/la_descriptor.h"
#    include "../tmp_config/met_flow_vars.h"
#    include "adaptivity/goal_functional.h"

///
/// \file goal_functional.h Abstract structure for Goal-functionals.
/// \author Martin Baumann, Philipp Gerstner
///
/// For a goal-functional \f$J\f$ of the following form
///   \f[
///        J(u):=\int_0^T (j_1(t,u(t)),\varphi)_\Omega dt + 
///                       (j_2(T,u(T)),\varphi)_\Omega.
///   \f]
/// The part \f$j_1\f$ is denoted by the function j_force_type and its related
/// parameters are stored in the struct force_type_parameters.
///

using hiflow::Vec;
using namespace hiflow::la;

namespace hiflow
{

    // **************************************************************************
    // Heat transfer in radial direction
    // **************************************************************************

    template<int DIM, class DataType>
    class GoalFunctionalHeatTransfer : public GoalFunctional<DIM, DataType>
    {
      public:
        GoalFunctionalHeatTransfer ( );
        void set_surface_type ( int type, int surface );
        DataType j_force_type ( ParametersForceType<DIM, DataType> p ) const;
        DataType j_final_type ( ParametersFinalType<DIM, DataType> p ) const;

        DataType j_force_eval ( ParametersEvalType<DIM, DataType> p ) const;
        DataType j_final_eval ( ParametersEvalType<DIM, DataType> p ) const;

        DataType j_force_eval_deriv ( ParametersEvalType<DIM, DataType> p ) const
        {
            return 0.;
        };

        DataType j_final_eval_deriv ( ParametersEvalType<DIM, DataType> p ) const
        {
            return 0.;
        };

      protected:
        int surface_;

    };

    // **************************************************************************
    // Kinetic Energy over regions where energy is greater than bound.
    // **************************************************************************

    template<int DIM, class DataType>
    class GoalFunctionalKineticEnergy : public GoalFunctional<DIM, DataType>
    {
      public:
        GoalFunctionalKineticEnergy ( );
        DataType j_force_type ( ParametersForceType<DIM, DataType> p ) const;
        DataType j_final_type ( ParametersFinalType<DIM, DataType> p ) const;

        DataType j_force_eval ( ParametersEvalType<DIM, DataType> p ) const;
        DataType j_final_eval ( ParametersEvalType<DIM, DataType> p ) const;

        DataType j_force_eval_deriv ( ParametersEvalType<DIM, DataType> p ) const
        {
            return 0.;
        };

        DataType j_final_eval_deriv ( ParametersEvalType<DIM, DataType> p ) const
        {
            return 0.;
        };

        DataType& bound ( )
        {
            return bound_;
        }

        void set_bound ( DataType low )
        {
            this->bound_ = low;
        }
      protected:
        DataType bound_;
    };

    // **************************************************************************
    // Vorticity over regions where vorticty is between bounds.
    // **************************************************************************

    template<int DIM, class DataType>
    class GoalFunctionalVorticity : public GoalFunctional<DIM, DataType>
    {
      public:
        GoalFunctionalVorticity ( );
        DataType j_force_type ( ParametersForceType<DIM, DataType> p ) const;
        DataType j_final_type ( ParametersFinalType<DIM, DataType> p ) const;

        DataType j_force_eval ( ParametersEvalType<DIM, DataType> p ) const;
        DataType j_final_eval ( ParametersEvalType<DIM, DataType> p ) const;

        DataType j_force_eval_deriv ( ParametersEvalType<DIM, DataType> p ) const
        {
            return 0.;
        };

        DataType j_final_eval_deriv ( ParametersEvalType<DIM, DataType> p ) const
        {
            return 0.;
        };

        void set_bounds ( DataType low, DataType high )
        {
            this->upper_bound_ = high;
            this->lower_bound_ = low;
        }

        void set_component ( int comp )
        {
            this->comp_ = comp;
        }

        void set_squared ( int squared )
        {
            this->squared_ = squared;
        }

        DataType& upper_bound ( )
        {
            return upper_bound_;
        }

        DataType& lower_bound ( )
        {
            return lower_bound_;
        }

      protected:
        int comp_;
        int squared_;
        ;
        DataType upper_bound_;
        DataType lower_bound_;
    };
}

#endif
