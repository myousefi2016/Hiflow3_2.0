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

#include "goal_functional.h"

#include <iostream>
#include <cmath>
#include <limits>

namespace hiflow
{

    template<int DIM, class DataType>
    GoalFunctional<DIM, DataType>::GoalFunctional ( )
    {
        force_type_active_ = false;
        final_type_active_ = false;
        name_ = "NOT SET";
        this->scale_ = 1.;
        this->co_system_ = 0;
        this->t_min_ = 0.;
        this->t_max_ = 99999.;
    }

    template<int DIM, class DataType>
    bool GoalFunctional<DIM, DataType>::is_in_box ( Vec<DIM, DataType> x ) const
    {
        switch ( this->co_system_ )
        {
            case 0:
                if ( x[0] < this->x_min_ || x[0] > this->x_max_ )
                    return false;

                if ( DIM > 1 )
                {
                    if ( x[1] < this->y_min_ || x[1] > this->y_max_ )
                        return false;

                    if ( DIM > 2 )
                    {
                        if ( x[2] < this->z_min_ || x[2] > this->z_max_ )
                            return false;
                    }
                }
                break;
            case 1:
                if ( x[0] < this->phi_min_ || x[0] > this->phi_max_ )
                    return false;

                if ( x[1] < this->r_min_ || x[1] > this->r_max_ )
                    return false;

                if ( DIM > 2 )
                {
                    if ( x[2] < this->z_min_ || x[2] > this->z_max_ )
                        return 0.;
                }
                break;
        }
        return true;
    }

    template class GoalFunctional<2, double>;
    template class GoalFunctional<3, double>;
    // **************************************************************************
    // Certain Variable in specifc subdomain of cylindrical geometry inside specific time interval
    // **************************************************************************

    template<int DIM, class DataType>
    GoalFunctionalVariableSubDomain<DIM, DataType>::GoalFunctionalVariableSubDomain ( )
    {
        this->name_ = "VariableSubDomain";

        this->force_type_active_ = true;
        this->final_type_active_ = true;
        this->integral_type_ = 1;
    }

    template<int DIM, class DataType>
    DataType GoalFunctionalVariableSubDomain<DIM, DataType>::j_force_type ( ParametersForceType<DIM, DataType> p ) const
    {
        if ( !this->force_type_active_ )
            return 0.;
        if ( !this->is_in_box ( p.x ) )
            return 0.;
        if ( p.absolute_time < this->t_min_ || p.absolute_time > this->t_max_ )
            return 0.;
        if ( p.var != var_ )
            return 0.;

        return p.phi * this->scale_;
    }

    template<int DIM, class DataType>
    DataType GoalFunctionalVariableSubDomain<DIM, DataType>::j_force_eval ( ParametersEvalType<DIM, DataType> p ) const
    {
        if ( !this->force_type_active_ )
            return 0.;
        if ( !this->is_in_box ( p.x ) )
            return 0.;
        if ( p.absolute_time < this->t_min_ || p.absolute_time > this->t_max_ )
            return 0.;

        return p.solP[var_] * this->scale_;
    }

    template<int DIM, class DataType>
    DataType GoalFunctionalVariableSubDomain<DIM, DataType>::j_final_type ( ParametersFinalType<DIM, DataType> p ) const
    {
        if ( !this->final_type_active_ )
            return 0.;
        if ( !this->is_in_box ( p.x ) )
            return 0.;
        if ( p.var != var_ )
            return 0.;

        return p.phi * this->scale_;
    }

    template<int DIM, class DataType>
    DataType GoalFunctionalVariableSubDomain<DIM, DataType>::j_final_eval ( ParametersEvalType<DIM, DataType> p ) const
    {
        if ( !this->final_type_active_ )
            return 0.;
        if ( !this->is_in_box ( p.x ) )
            return 0.;

        return p.solP[var_] * this->scale_;
    }

    template<int DIM, class DataType>
    DataType GoalFunctionalVariableSubDomain<DIM, DataType>::j_force_eval_deriv ( ParametersEvalType<DIM, DataType> p ) const
    {
        if ( !this->force_type_active_ )
            return 0.;
        if ( !this->is_in_box ( p.x ) )
            return 0.;
        if ( p.absolute_time < this->t_min_ || p.absolute_time > this->t_max_ )
            return 0.;
        if ( p.var != var_ )
            return 0.;

        return p.d_solP[var_] * this->scale_;
    }

    template<int DIM, class DataType>
    DataType GoalFunctionalVariableSubDomain<DIM, DataType>::j_final_eval_deriv ( ParametersEvalType<DIM, DataType> p ) const
    {
        if ( !this->final_type_active_ )
            return 0.;
        if ( !this->is_in_box ( p.x ) )
            return 0.;
        if ( p.var != var_ )
            return 0.;

        return p.d_solP[var_] * this->scale_;
    }

    template class GoalFunctionalVariableSubDomain<2, double>;
    template class GoalFunctionalVariableSubDomain<3, double>;

}
