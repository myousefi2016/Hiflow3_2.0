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

#ifndef HIFLOW_ADAPTIVITY_GOAL_FUNCTIONAL_H
#    define HIFLOW_ADAPTIVITY_GOAL_FUNCTIONAL_H

#    include <iostream>
#    include <string>
#    include "assembly/function_values.h"
#    include "common/vector_algebra.h"
#    include "assembly/function_values.h"
#    include "linear_algebra/la_descriptor.h"

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

    template<int DIM, class DataType>
    struct ParametersForceType
    {
        Vec<DIM, DataType> x;
        DataType absolute_time;
        int var;
        DataType relative_time; // -1 -> t=t_{i-1}, 0 -> t=t_i, 1->t=t_{i+1}
        std::vector<DataType> solP_prev; // primal solution t_{i-1}
        std::vector<DataType> solP; // primal solution t_i
        std::vector<DataType> solP_next; // primal solution t_i+1
        std::vector< Vec<DIM, DataType> > grad_solP_prev; // primal solution t_{i-1}
        std::vector< Vec<DIM, DataType> > grad_solP; // primal solution t_i
        std::vector< Vec<DIM, DataType> > grad_solP_next; // primal solution t_i+1
        DataType phi; // test function for variable 'var'
        Vec<DIM, DataType> grad_phi; // gradient of test function for variable 'var'
    };

    template<int DIM, class DataType>
    struct ParametersFinalType
    {
        Vec<DIM, DataType> x;
        int var;
        std::vector<DataType> solP; // primal solution
        std::vector< Vec<DIM, DataType> > grad_solP; // primal solution
        DataType phi; // test function for variable 'var'
        Vec<DIM, DataType> grad_phi; // gradient of test function for variable 'var'
    };

    template<int DIM, class DataType>
    struct ParametersEvalType
    {
        Vec<DIM, DataType> x;
        DataType absolute_time;
        int var;
        std::vector<DataType> solP; // primal solution t_i
        std::vector< Vec<DIM, DataType> > grad_solP; // primal solution t_i
        std::vector<DataType> d_solP; // primal solution t_i
        std::vector< Vec<DIM, DataType> > grad_d_solD; // primal solution t_i
    };

    template<int DIM, class DataType >
    class GoalFunctional
    {
      public:

        GoalFunctional ( );

        virtual ~GoalFunctional ( )
        {
            ;
        }

        /// \brief Force-type contribution of a goal-functional
        /// Only the part \f$(j_1(t,u(t)),varphi)\f$ should be implemented, excluding
        /// parts related to quadrature in space or in time.
        virtual DataType j_force_type ( ParametersForceType<DIM, DataType> p ) const = 0;

        /// \brief final-type contribution of a goal-functional (related to \f$t=T\f$)
        /// Only the part \f$(j_1(t,u(t)),varphi)\f$ should be implemented, excluding
        /// parts related to quadrature in space or in time.
        virtual DataType j_final_type ( ParametersFinalType<DIM, DataType> p ) const = 0;

        virtual DataType j_force_eval ( ParametersEvalType<DIM, DataType> p ) const = 0;

        virtual DataType j_final_eval ( ParametersEvalType<DIM, DataType> p ) const = 0;

        virtual DataType j_force_eval_deriv ( ParametersEvalType<DIM, DataType> p ) const
        {
            return 0.;
        };

        virtual DataType j_final_eval_deriv ( ParametersEvalType<DIM, DataType> p ) const
        {
            return 0.;
        };

        bool& force_type_active ( )
        {
            return force_type_active_;
        }

        bool& final_type_active ( )
        {
            return final_type_active_;
        }

        void set_active_parts ( bool active_force, bool active_ic )
        {
            this->force_type_active_ = active_force;
            this->final_type_active_ = active_ic;
        }

        void set_cyl_box ( DataType t_min, DataType t_max, DataType z_min, DataType z_max, DataType r_min, DataType r_max, DataType phi_min, DataType phi_max )
        {
            this->t_min_ = t_min;
            this->t_max_ = t_max;
            this->z_min_ = z_min;
            this->z_max_ = z_max;
            this->r_min_ = r_min;
            this->r_max_ = r_max;
            this->phi_min_ = phi_min;
            this->phi_max_ = phi_max;
            this->co_system_ = 1;
        }

        void set_cart_box ( DataType t_min, DataType t_max, DataType x_min, DataType x_max, DataType y_min, DataType y_max, DataType z_min, DataType z_max )
        {
            this->t_min_ = t_min;
            this->t_max_ = t_max;
            this->z_min_ = z_min;
            this->z_max_ = z_max;
            this->y_min_ = y_min;
            this->y_max_ = y_max;
            this->x_min_ = x_min;
            this->x_max_ = x_max;
            this->co_system_ = 0;
        }

        void set_scale_factor ( DataType scale )
        {
            this->scale_ = scale;
        }

        int GetIntegralType ( )
        {
            return this->integral_type_;
        }

        void use_cart_co_system ( )
        {
            this->co_system_ = 0;
        }

        void use_cyl_co_system ( )
        {
            this->co_system_ = 1;
        }

      protected:

        bool is_in_box ( Vec<DIM, DataType> x ) const;

        std::string name_;
        bool force_type_active_;
        bool final_type_active_;

        DataType x_min_;
        DataType x_max_;
        DataType y_min_;
        DataType y_max_;
        DataType z_min_;
        DataType z_max_;
        DataType r_min_;
        DataType r_max_;
        DataType phi_min_;
        DataType phi_max_;
        DataType t_min_;
        DataType t_max_;
        DataType scale_;

        int co_system_;
        int integral_type_;
    };

    // **************************************************************************
    // Specific veriable in specifc subdomain of cylindrical geometry during specific time interval
    // **************************************************************************

    template<int DIM, class DataType>
    class GoalFunctionalVariableSubDomain : public GoalFunctional<DIM, DataType>
    {
      public:
        GoalFunctionalVariableSubDomain ( );

        void set_variable ( int var )
        {
            this->var_ = var;
        }

        DataType j_force_type ( ParametersForceType<DIM, DataType> p ) const;
        DataType j_final_type ( ParametersFinalType<DIM, DataType> p ) const;

        DataType j_force_eval ( ParametersEvalType<DIM, DataType> p ) const;
        DataType j_final_eval ( ParametersEvalType<DIM, DataType> p ) const;
        DataType j_force_eval_deriv ( ParametersEvalType<DIM, DataType> p ) const;
        DataType j_final_eval_deriv ( ParametersEvalType<DIM, DataType> p ) const;

      protected:

        int var_;
    };
}

#endif
