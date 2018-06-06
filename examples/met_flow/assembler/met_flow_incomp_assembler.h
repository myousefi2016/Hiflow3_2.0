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

#ifndef MET_FLOW_INCOMP_ASSEMBLER_H
#    define MET_FLOW_INCOMP_ASSEMBLER_H

///
/// \brief Assembler base class for incompressible Navier Stokes equations
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
#    include "../general/goal_functional.h"
#    include "../met_flow_assembler.h"

template<int DIM, class DataType>
class MetFlowIncompAssembler : public virtual MetFlowAssembler<DIM, DataType>
{
  public:

    MetFlowIncompAssembler ( );

    ~MetFlowIncompAssembler ( )
    {
        ;
    }

    //////////////////////////////////////////
    //// Functions from MetFlowAssembler /////
    virtual void set_time_discretization_to_zero ( );
    virtual void set_time_discretization_to_simple ( );
    virtual void set_time_discretization_to_theta ( DataType theta );
    virtual void set_time_discretization_to_galerkin_cd ( bool modified );
    virtual void set_time_discretization_to_galerkin_dc ( bool modified );
    virtual void set_time_discretization_off ( );
    virtual void set_time_discretization_to_scaling ( DataType beta1, DataType beta2, DataType beta3, DataType beta4, DataType beta5, DataType beta6 );

    virtual void print_parameters ( );
    virtual void print_coeff ( );
    virtual DataType get_coeff ( std::string term );

    //////////////////////////////////////////////////////
    ///// Specific functions for incompressible flow /////

    void set_graddiv ( DataType gamma, DataType C )
    {
        graddiv_gamma_ = gamma;
        graddiv_C_ = C;
    }

    void set_omega ( DataType omega )
    {
        omega_ = omega;
    }

    void set_nu ( DataType nu )
    {
        nu_ = nu;
    }

    void set_rho ( DataType rho )
    {
        interminable_assert ( fabs ( rho ) > 1.0e-15 );
        inv_rho_ = 1. / rho;
        rho_ = rho;
    }

    void set_conv_mode_to_stokes ( )
    {
        this->conv_mode_ = STOKES;
    }

    void set_conv_mode_to_oseen ( )
    {
        this->conv_mode_ = OSEEN;
    }

    void set_conv_mode_to_navier_stokes ( )
    {
        this->conv_mode_ = NAVIERSTOKES;
    }

    void set_vector_mode_to_std ( )
    {
        this->vector_asm_mode_ = VECTOR_STD;
    }

    void set_vector_mode_to_goal ( )
    {
        this->vector_asm_mode_ = VECTOR_GOAL;
    }

    void set_vector_mode_to_vort ( )
    {
        this->vector_asm_mode_ = VECTOR_VORT;
    }

    void set_vector_mode_to_quant ( int incomp_scalars )
    {
        this->vector_asm_mode_ = VECTOR_QUANT;
        this->incomp_scalars_ = incomp_scalars;
        this->num_scalars_ = incomp_scalars;
    }

    void set_scalar_mode_to_div_mean ( )
    {
        this->sca_mode_ = DIV_MEAN;
    }

    void set_scalar_mode_to_div_max ( )
    {
        this->sca_mode_ = DIV_MAX;
    }

    void set_graddiv_mode_to_OFF ( )
    {
        graddiv_mode_ = GD_OFF;
    }

    void set_graddiv_mode_to_CONST ( )
    {
        graddiv_mode_ = CONST;
    }

    void set_graddiv_mode_to_VMS ( )
    {
        graddiv_mode_ = VMS;
    }

    void set_graddiv_mode_to_SUPG ( )
    {
        graddiv_mode_ = SUPG;
    }

    void set_skew_mode_to_OFF ( )
    {
        skew_mode_ = SS_OFF;
    }

    void set_skew_mode_to_ON ( )
    {
        skew_mode_ = ON;
    }

    virtual void clear ( );

  protected:

    ///// Specific incompressible flow functions /////

    virtual void assemble_local_matrix_primal ( const Element<DataType>& element, LocalMatrix& lm ) const
    {
        ;
    }

    virtual void assemble_local_matrix_dual ( const Element<DataType>& element, LocalMatrix& lm ) const
    {
        ;
    }

    virtual void assemble_local_vector_primal ( const Element<DataType>& element, LocalVector& lv ) const
    {
        ;
    }

    virtual void assemble_local_vector_dual ( const Element<DataType>& element, LocalVector& lv ) const
    {
        ;
    }

    virtual void assemble_local_vector_vort ( const Element<DataType>& element, LocalVector& lv ) const
    {
        ;
    }

    virtual void assemble_local_vector_quant ( const Element<DataType>& element, LocalVector& lv ) const
    {
        ;
    }

    virtual void assemble_local_vector_goal ( const Element<DataType>& element, LocalVector& lv ) const
    {
        ;
    }

    virtual void assemble_local_scalar_goal_int ( const Element<DataType>& element, DataType& ls ) const
    {
        ;
    }

    virtual void assemble_local_scalar_goal_fin ( const Element<DataType>& element, DataType& ls ) const
    {
        ;
    }

    virtual void assemble_local_scalar_div_max ( const Element<DataType>& element, DataType& ls ) const
    {
        ;
    }

    virtual void assemble_local_scalar_div_mean ( const Element<DataType>& element, DataType& ls ) const
    {
        ;
    }

    virtual void compute_tau_graddiv ( const Element<DataType>& element, DataType& tau ) const;

    // fluid parameters
    DataType nu_, inv_rho_, rho_;
    DataType omega_;

    // different modes
    GradDivMode graddiv_mode_;
    SkewMode skew_mode_;
    VectorMode vector_asm_mode_;
    ConvectionMode conv_mode_;

    DataType graddiv_gamma_;
    DataType graddiv_C_;

    int incomp_scalars_;
    int num_scalars_;

    // time stepping coefficients
    DataType theta_mom_vis_c_; // nu * Laplace(u_n)
    DataType theta_mom_vis_p_; // nu * Laplace(u_(n-1))
    DataType theta_mom_adv_cc_; // u_n * nabla(u_n)
    DataType theta_mom_adv_pc_; // u_(n-1) * nabla(u_n)
    DataType theta_mom_adv_cp_; // u_n * nabla(u_(n-1))
    DataType theta_mom_adv_pp_; // u_(n-1) * nabla(u_(n-1))
    DataType theta_mom_pre_c_; // grad(p_n)
    DataType theta_mom_graddiv_c_; // div(u_n)*div(v)
    DataType theta_mom_graddiv_p_; // div(u_(n-1))*div(v)
    DataType theta_mom_rot_c_; // Cor(u_n)
    DataType theta_mom_rot_p_; // Cor(u_(n-1))
    DataType theta_inc_c_; // div(u_n)
    DataType theta_inc_p_; // div(u_(n-1))

    DataType delta_mass_;
    DataType delta_mom_vis_c_; // nu * Laplace(u_n)
    DataType delta_mom_vis_n_; // nu * Laplace(u_(n-1))
    DataType delta_mom_adv_cc_; // u_n * nabla(u_n)
    DataType delta_mom_adv_nc_; // u_(n-1) * nabla(u_n)
    DataType delta_mom_adv_cn_; // u_n * nabla(u_(n-1))
    DataType delta_mom_adv_nn_; // u_(n-1) * nabla(u_(n-1))
    DataType delta_mom_adv_pc_; // u_(n-1) * nabla(u_(n-1))
    DataType delta_mom_pre_c_; // grad(p_n)
    DataType delta_mom_pre_n_; // grad(p_n)
    DataType delta_mom_graddiv_c_; // div(u_n)*div(v)
    DataType delta_mom_graddiv_n_; // div(u_(n-1))*div(v)
    DataType delta_mom_rot_c_; // Cor(u_n)
    DataType delta_mom_rot_n_; // Cor(u_(n-1))
    DataType delta_mom_rea_cc_;
    DataType delta_mom_rea_nc_;
    DataType delta_mom_rea_cn_;
    DataType delta_mom_rea_nn_;
    DataType delta_mom_rea_pc_;

    DataType delta_inc_c_; // div(u_n)
    DataType delta_inc_n_; // div(u_(n-1))

    using MetFlowAssembler<DIM, DataType>::num_dofs;
    using MetFlowAssembler<DIM, DataType>::dof_index;
    using MetFlowAssembler<DIM, DataType>::phi;
    using MetFlowAssembler<DIM, DataType>::grad_phi;
};

#endif
