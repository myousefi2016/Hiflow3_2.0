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

    // **************************************************************************
    // Heat transfer in radial direction through the [rMin,rMax]
    // **************************************************************************

    template<int DIM, class DataType>
    GoalFunctionalHeatTransfer<DIM, DataType>::GoalFunctionalHeatTransfer ( )
    {
        this->name_ = "HeatTransfer";

        this->force_type_active_ = true;
        this->final_type_active_ = true;
        this->integral_type_ = 1;
    }

    template<int DIM, class DataType>
    void GoalFunctionalHeatTransfer<DIM, DataType>::set_surface_type ( int type, int surface )
    {
        this->integral_type_ = type;
        this->surface_ = surface;
    }

    template<int DIM, class DataType>
    DataType GoalFunctionalHeatTransfer<DIM, DataType>::j_force_type ( ParametersForceType<DIM, DataType> p ) const
    {
        if ( !this->force_type_active_ )
            return 0.;
        if ( !this->is_in_box ( p.x ) )
            return 0.;

        if ( p.absolute_time < this->t_min_ || p.absolute_time > this->t_max_ )
            return 0.;

        DataType ret = 0.;
        if ( this->integral_type_ == 1 )
        {
            if ( p.var == DIM + 1 )
            {
                ret += p.grad_phi[1];
            }
        }
        else
        {

        }
        return ret * this->scale_;
    }

    template<int DIM, class DataType>
    DataType GoalFunctionalHeatTransfer<DIM, DataType>::j_force_eval ( ParametersEvalType<DIM, DataType> p ) const
    {
        if ( !this->force_type_active_ )
            return 0.;
        if ( !this->is_in_box ( p.x ) )
            return 0.;

        if ( p.absolute_time < this->t_min_ || p.absolute_time > this->t_max_ )
            return 0.;

        DataType ret = 0.;
        if ( this->integral_type_ == 1 )
        {
            ret += p.grad_solP[DIM + 1][1];
        }
        else
        {

        }
        return ret * this->scale_;
    }

    template<int DIM, class DataType>
    DataType GoalFunctionalHeatTransfer<DIM, DataType>::j_final_type ( ParametersFinalType<DIM, DataType> p ) const
    {
        if ( !this->final_type_active_ )
            return 0.;
        if ( !this->is_in_box ( p.x ) )
            return 0.;

        DataType ret = 0.;
        if ( this->integral_type_ == 1 )
        {
            if ( p.var == DIM + 1 )
            {
                ret += p.grad_phi[1];
            }
        }
        else
        {

        }

        return ret * this->scale_;
    }

    template<int DIM, class DataType>
    DataType GoalFunctionalHeatTransfer<DIM, DataType>::j_final_eval ( ParametersEvalType<DIM, DataType> p ) const
    {
        if ( !this->final_type_active_ )
            return 0.;
        if ( !this->is_in_box ( p.x ) )
            return 0.;

        DataType ret = 0.;
        if ( this->integral_type_ == 1 )
        {
            ret += p.grad_solP[DIM + 1][1];
        }
        else
        {

        }

        return ret * this->scale_;
    }

    template class GoalFunctionalHeatTransfer<2, double>;
    template class GoalFunctionalHeatTransfer<3, double>;

    // **************************************************************************
    // Kinetic Energy over regions where energy is greater than bound.
    // **************************************************************************

    /// \brief Kinetic integrated over regions where energy is greater than bound.
    /// Definition of Goal-functional:
    ///   \f[ J(u) := \int_0^T \frac{1}{2} (\Chi(t)\cdot u(t),u(t))_\Omega dt + 
    ///                        \frac{1}{2} (\Chi(T)\cdot u(T),u(T))_\Omega. \f]
    /// Contribution to dual problem:
    ///   \f[ \nabla_uJ(u)\varphi := \int_0^T (\Chi(t)\cdot u(t),\varphi)_\Omega dt + 
    ///                                       (\Chi(T)\cdot u(T),\varphi)_\Omega. \f]

    template<int DIM, class DataType>
    GoalFunctionalKineticEnergy<DIM, DataType>::GoalFunctionalKineticEnergy ( )
    {
        this->name_ = "KineticEnergy";

        // bounds set to be global integral
        bound_ = 0.;
        this->integral_type_ = 1;
    }

    template<int DIM, class DataType>
    DataType GoalFunctionalKineticEnergy<DIM, DataType>::j_force_type ( ParametersForceType<DIM, DataType> p ) const
    {
        if ( !this->force_type_active_ )
            return 0.;
        if ( !this->is_in_box ( p.x ) )
            return 0.;
        if ( p.absolute_time < this->t_min_ || p.absolute_time > this->t_max_ )
            return 0.;
        if ( p.var >= DIM )
            return 0.;

        std::vector<DataType> solP;
        std::vector< Vec<DIM, DataType> > grad_solP;

        if ( p.relative_time == -1 )
        {
            solP = p.solP_prev;
            grad_solP = p.grad_solP_prev;
        }
        if ( p.relative_time == -0.5 )
        {
            solP = p.solP_prev;
            grad_solP = p.grad_solP_prev;
            for ( int l = 0; l < solP.size ( ); ++l )
            {
                solP[l] += p.solP[l];
                solP[l] *= 0.5;
                for ( int d = 0; d < DIM; ++d )
                {
                    grad_solP[l][d] += p.grad_solP[l][d];
                    grad_solP[l][d] *= 0.5;
                }
            }
        }
        if ( p.relative_time == 0 )
        {
            solP = p.solP;
            grad_solP = p.grad_solP;
        }
        if ( p.relative_time == 0.5 )
        {
            solP = p.solP_next;
            grad_solP = p.grad_solP_next;
            for ( int l = 0; l < solP.size ( ); ++l )
            {
                solP[l] += p.solP[l];
                solP[l] *= 0.5;
                for ( int d = 0; d < DIM; ++d )
                {
                    grad_solP[l][d] += p.grad_solP[l][d];
                    grad_solP[l][d] *= 0.5;
                }
            }
        }
        if ( p.relative_time == 1 )
        {
            solP = p.solP_next;
            grad_solP = p.grad_solP_next;
        }

        DataType ret = 0.;
        DataType energy = solP[0] * solP[0] + solP[1] * solP[1];
        if ( DIM == 3 )
            energy += solP[2] * solP[2];

        if ( energy > bound_ )
            ret += 2.0 * solP[p.var] * p.phi;

        return ret * this->scale_;
    }

    template<int DIM, class DataType>
    DataType GoalFunctionalKineticEnergy<DIM, DataType>::j_force_eval ( ParametersEvalType<DIM, DataType> p ) const
    {
        if ( !this->force_type_active_ )
            return 0.;
        if ( !this->is_in_box ( p.x ) )
            return 0.;
        if ( p.absolute_time < this->t_min_ || p.absolute_time > this->t_max_ )
            return 0.;
        if ( p.var >= DIM )
            return 0.;

        DataType ret = 0.;
        DataType energy = p.solP[0] * p.solP[0] + p.solP[1] * p.solP[1];
        if ( DIM == 3 )
            energy += p.solP[2] * p.solP[2];

        if ( energy > bound_ )
            ret += energy;

        return ret * this->scale_;
    }

    template<int DIM, class DataType>
    DataType GoalFunctionalKineticEnergy<DIM, DataType>::j_final_type ( ParametersFinalType<DIM, DataType> p ) const
    {
        if ( !this->final_type_active_ )
            return 0.;
        if ( !this->is_in_box ( p.x ) )
            return 0.;
        if ( p.var >= DIM )
            return 0.;

        DataType energy = ( p.solP )[0]*( p.solP )[0] + ( p.solP )[1]*( p.solP )[1];
        if ( DIM == 3 )
            energy += ( p.solP )[2]*( p.solP )[2];

        DataType ret = 0.;
        if ( energy > bound_ )
            ret += 2.0 * p.solP[p.var] * p.phi;

        return ret * this->scale_;
    }

    template<int DIM, class DataType>
    DataType GoalFunctionalKineticEnergy<DIM, DataType>::j_final_eval ( ParametersEvalType<DIM, DataType> p ) const
    {
        if ( !this->final_type_active_ )
            return 0.;
        if ( !this->is_in_box ( p.x ) )
            return 0.;
        if ( p.var >= DIM )
            return 0.;

        DataType energy = ( p.solP )[0]*( p.solP )[0] + ( p.solP )[1]*( p.solP )[1];
        if ( DIM == 3 )
            energy += ( p.solP )[2]*( p.solP )[2];

        DataType ret = 0.;
        if ( energy > bound_ )
            ret += energy;

        return ret * this->scale_;
    }

    template class GoalFunctionalKineticEnergy<2, double>;
    template class GoalFunctionalKineticEnergy<3, double>;

    // **************************************************************************
    // Vorticity over regions where vorticty is between bounds.
    // **************************************************************************

    /// \brief Vorticity integrated over regions where vorticity is between bounds.
    /// Vorticity corresponds to the scalar valued horizontal vorticity in case
    /// of 3D, i.e. \f$\zeta:=\partial_x v-\partial_y u\f$. 
    /// Definition of Goal-functional:
    ///   \f[ J(u) := \int_0^T \frac{1}{2} (\Chi(t),\nabla \times u(t))_\Omega dt + 
    ///                        \frac{1}{2} (\Chi(T),\nabla \times u(T))_\Omega. \f]
    /// Contribution to dual problem:
    ///   \f[ \nabla_uJ(u)\varphi := \int_0^T (\Chi(t),\nabla \times \varphi)_\Omega dt + 
    ///                                       (\Chi(T),\nabla \times \varphi)_\Omega. \f]

    template<int DIM, class DataType>
    GoalFunctionalVorticity<DIM, DataType>::GoalFunctionalVorticity ( )
    {
        this->name_ = "Vorticity";

        // bounds set to be global integral
        lower_bound_ = std::numeric_limits<DataType>::min ( );
        upper_bound_ = std::numeric_limits<DataType>::max ( );
        this->integral_type_ = 1;
    }

    template<int DIM, class DataType>
    DataType GoalFunctionalVorticity<DIM, DataType>::j_force_type ( ParametersForceType<DIM, DataType> p ) const
    {
        if ( !this->force_type_active_ )
            return 0.;
        if ( !this->is_in_box ( p.x ) )
            return 0.;
        if ( p.absolute_time < this->t_min_ || p.absolute_time > this->t_max_ )
            return 0.;

        std::vector<DataType> solP;
        std::vector< Vec<DIM, DataType> > grad_solP;

        if ( p.relative_time == -1 )
        {
            solP = p.solP_prev;
            grad_solP = p.grad_solP_prev;
        }
        if ( p.relative_time == -0.5 )
        {
            solP = p.solP_prev;
            grad_solP = p.grad_solP_prev;
            for ( int l = 0; l < solP.size ( ); ++l )
            {
                solP[l] += p.solP[l];
                solP[l] *= 0.5;
                for ( int d = 0; d < DIM; ++d )
                {
                    grad_solP[l][d] += p.grad_solP[l][d];
                    grad_solP[l][d] *= 0.5;
                }
            }
        }
        if ( p.relative_time == 0 )
        {
            solP = p.solP;
            grad_solP = p.grad_solP;
        }
        if ( p.relative_time == 0.5 )
        {
            solP = p.solP_next;
            grad_solP = p.grad_solP_next;
            for ( int l = 0; l < solP.size ( ); ++l )
            {
                solP[l] += p.solP[l];
                solP[l] *= 0.5;
                for ( int d = 0; d < DIM; ++d )
                {
                    grad_solP[l][d] += p.grad_solP[l][d];
                    grad_solP[l][d] *= 0.5;
                }
            }
        }
        if ( p.relative_time == 1 )
        {
            solP = p.solP_next;
            grad_solP = p.grad_solP_next;
        }

        DataType ret = 0.;
        DataType inv_r = 1. / p.x[1];

        DataType vorticity = 0.;

        switch ( this->comp_ )
        {
            case 0:
                vorticity = grad_solP[1][2] - grad_solP[2][1];
                break;
            case 1:
                vorticity = inv_r * grad_solP[2][0] - grad_solP[0][2];
                break;
            case 2:
                vorticity = inv_r * ( solP[0] + p.x[1] * grad_solP[0][1] - grad_solP[1][0] );
                break;
        }
        DataType abs_vorticity;
        if ( this->squared_ == 1 )
            abs_vorticity = vorticity * vorticity;
        else
            abs_vorticity = vorticity;

        if ( ( abs_vorticity >= lower_bound_ ) && ( abs_vorticity <= upper_bound_ ) )
        {
            if ( p.var == 0 )
            {
                if ( this->comp_ == 1 )
                    ret = -p.grad_phi[2];
                if ( this->comp_ == 2 )
                    ret = inv_r * ( p.phi + p.x[1] * p.grad_phi[1] );
            }
            if ( p.var == 1 )
            {
                if ( this->comp_ == 0 )
                    ret = p.grad_phi[2];
                if ( this->comp_ == 2 )
                    ret = -inv_r * p.grad_phi[0];
            }
            if ( p.var == 2 )
            {
                if ( this->comp_ == 0 )
                    ret = -p.grad_phi[1];
                if ( this->comp_ == 1 )
                    ret = inv_r * p.grad_phi[0];
            }
        }
        if ( this->squared_ == 1 )
        {
            ret *= 2. * vorticity;
        }

        return ret * this->scale_;
    }

    template<int DIM, class DataType>
    DataType GoalFunctionalVorticity<DIM, DataType>::j_force_eval ( ParametersEvalType<DIM, DataType> p ) const
    {
        if ( !this->force_type_active_ )
            return 0.;
        if ( !this->is_in_box ( p.x ) )
            return 0.;
        if ( p.absolute_time < this->t_min_ || p.absolute_time > this->t_max_ )
            return 0.;

        DataType ret = 0.;
        DataType inv_r = 1. / p.x[1];

        DataType vorticity = 0.;

        switch ( this->comp_ )
        {
            case 0:
                vorticity = p.grad_solP[1][2] - p.grad_solP[2][1];
                break;
            case 1:
                vorticity = inv_r * p.grad_solP[2][0] - p.grad_solP[0][2];
                break;
            case 2:
                vorticity = inv_r * ( p.solP[0] + p.x[1] * p.grad_solP[0][1] - p.grad_solP[1][0] );
                break;
        }
        DataType abs_vorticity;
        if ( this->squared_ == 1 )
            abs_vorticity = vorticity * vorticity;
        else
            abs_vorticity = vorticity;

        if ( ( abs_vorticity >= lower_bound_ ) && ( abs_vorticity <= upper_bound_ ) )
        {
            ret += abs_vorticity;
        }

        return ret * this->scale_;
    }

    template<int DIM, class DataType>
    DataType GoalFunctionalVorticity<DIM, DataType>::j_final_type ( ParametersFinalType<DIM, DataType> p ) const
    {
        if ( !this->final_type_active_ )
            return 0.;
        if ( !this->is_in_box ( p.x ) )
            return 0.;

        DataType ret = 0.;
        DataType inv_r = 1. / p.x[1];
        DataType vorticity = 0.;

        switch ( this->comp_ )
        {
            case 0:
                vorticity = p.grad_solP[1][2] - p.grad_solP[2][1];
                break;
            case 1:
                vorticity = inv_r * p.grad_solP[2][0] - p.grad_solP[0][2];
                break;
            case 2:
                vorticity = inv_r * ( p.solP[0] + p.x[1] * p.grad_solP[0][1] - p.grad_solP[1][0] );
                break;
        }
        DataType abs_vorticity;
        if ( this->squared_ == 1 )
            abs_vorticity = vorticity * vorticity;
        else
            abs_vorticity = vorticity;

        if ( ( abs_vorticity >= lower_bound_ ) && ( abs_vorticity <= upper_bound_ ) )
        {
            if ( p.var == 0 )
            {
                if ( this->comp_ == 1 )
                    ret = -p.grad_phi[2];
                if ( this->comp_ == 2 )
                    ret = inv_r * ( p.phi + p.x[1] * p.grad_phi[1] );
            }
            if ( p.var == 1 )
            {
                if ( this->comp_ == 0 )
                    ret = p.grad_phi[2];
                if ( this->comp_ == 2 )
                    ret = -inv_r * p.grad_phi[0];
            }
            if ( p.var == 2 )
            {
                if ( this->comp_ == 0 )
                    ret = -p.grad_phi[1];
                if ( this->comp_ == 1 )
                    ret = inv_r * p.grad_phi[0];
            }
        }

        if ( this->squared_ == 1 )
            ret *= 2. * vorticity;

        return ret * this->scale_;
    }

    template<int DIM, class DataType>
    DataType GoalFunctionalVorticity<DIM, DataType>::j_final_eval ( ParametersEvalType<DIM, DataType> p ) const
    {
        if ( !this->final_type_active_ )
            return 0.;
        if ( !this->is_in_box ( p.x ) )
            return 0.;

        DataType ret = 0.;
        DataType inv_r = 1. / p.x[1];
        DataType vorticity = 0.;

        switch ( this->comp_ )
        {
            case 0:
                vorticity = p.grad_solP[1][2] - p.grad_solP[2][1];
                break;
            case 1:
                vorticity = inv_r * p.grad_solP[2][0] - p.grad_solP[0][2];
                break;
            case 2:
                vorticity = inv_r * ( p.solP[0] + p.x[1] * p.grad_solP[0][1] - p.grad_solP[1][0] );
                break;
        }
        DataType abs_vorticity;
        if ( this->squared_ == 1 )
            abs_vorticity = vorticity * vorticity;
        else
            abs_vorticity = vorticity;

        if ( ( abs_vorticity >= lower_bound_ ) && ( abs_vorticity <= upper_bound_ ) )
        {
            ret += abs_vorticity;
        }

        return ret * this->scale_;
    }

    template class GoalFunctionalVorticity<2, double>;
    template class GoalFunctionalVorticity<3, double>;
}
