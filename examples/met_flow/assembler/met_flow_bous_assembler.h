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

#ifndef MET_FLOW_BOUS_ASSEMBLER_H
#    define MET_FLOW_BOUS_ASSEMBLER_H

///
/// \brief Assembler base class for Boussinesq equations
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
#    include "../met_flow_incomp_assembler.h"

template<int DIM, class DataType>
class MetFlowBousAssembler : public virtual MetFlowIncompAssembler<DIM, DataType>
{
  public:

    MetFlowBousAssembler ( );

    ~MetFlowBousAssembler ( )
    {
        ;
    }

    ////////////////////////////////////////////////
    //// Functions from MetFlowIncompAssembler /////
    virtual void set_time_discretization_to_simple ( );
    virtual void set_time_discretization_to_zero ( );
    virtual void set_time_discretization_to_theta ( DataType theta );
    virtual void set_time_discretization_to_galerkin_cd ( bool modified );
    virtual void set_time_discretization_to_galerkin_dc ( bool modified );
    virtual void set_time_discretization_off ( );

    virtual void print_parameters ( );
    virtual void print_coeff ( );
    virtual DataType get_coeff ( std::string term );
    virtual void initialize_for_element ( const Element<DataType>& element, const Quadrature<DataType>& quadrature );
    virtual void initialize_for_facet ( const Element<DataType>& element, const Quadrature<DataType>& quadrature, int facet_number );

    ////////////////////////////////////////////////
    //// Specific functions for Boussinesq /////////

    void set_kappa ( DataType d )
    {
        kappa_ = d;
    }

    void set_alpha_g ( DataType alpha )
    {
        alpha_g_ = alpha;
    }

    void set_ref_T ( DataType ref_T )
    {
        ref_T_ = ref_T;
    }

    void set_start_T ( DataType start_T )
    {
        start_T_ = start_T;
    }

    void set_temp_supg ( DataType gamma, DataType azi, DataType rad, DataType axi )
    {
        this->temp_supg_gamma_ = gamma;
        this->temp_supg_comp_.resize ( DIM, 0. );
        this->temp_supg_comp_[0] = azi;
        this->temp_supg_comp_[1] = rad;
        this->temp_supg_comp_[DIM - 1] = axi;
    }
    void set_gravity ( std::vector<DataType> g_cart );

    void set_nusselt_area ( DataType r_min, DataType r_max, int id )
    {
        nusselt_r_min_ = r_min;
        nusselt_r_max_ = r_max;
        nusselt_surface_id_ = id;
    }

    void set_vector_mode_to_grad_temp ( )
    {
        this->vector_asm_mode_ = VECTOR_GRADTEMP;
    }

    void set_scalar_mode_to_temp_mean ( )
    {
        this->sca_mode_ = TEMP_MEAN;
    }

    void set_scalar_mode_to_rad_heat_volume ( )
    {
        this->sca_mode_ = RAD_HEAT_VOLUME;
    }

    void set_scalar_mode_to_rad_heat_surface ( )
    {
        this->sca_mode_ = RAD_HEAT_SURFACE;
    }

    void set_temp_supg_mode_to_OFF ( )
    {
        this->temp_supg_mode_ = TS_OFF;
    }

    void set_temp_supg_mode_to_STD ( )
    {
        this->temp_supg_mode_ = STD;
    }

    void set_temp_supg_mode_to_DEP ( )
    {
        this->temp_supg_mode_ = DEP;
    }

    void set_temp_supg_mode_to_DEP_IMPLICIT ( )
    {
        this->temp_supg_mode_ = DEP_IMPLICIT;
    }

    void set_vector_mode_to_quant ( int incomp_scalars, int bous_scalars )
    {
        this->vector_asm_mode_ = VECTOR_QUANT;
        this->incomp_scalars_ = incomp_scalars;
        this->bous_scalars_ = bous_scalars;
        this->num_scalars_ = incomp_scalars + bous_scalars;
    }

    virtual void clear ( );

  protected:
    ////////////////////////////////////////////////
    //// Functions from MetFlowIncompAssembler /////

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

    virtual void assemble_local_vector_quant ( const Element<DataType>& element, LocalVector& lv ) const
    {
        ;
    }

    virtual void assemble_local_vector_goal ( const Element<DataType>& element, LocalVector& lv ) const
    {
        ;
    }

    virtual void assemble_local_vector_buoyancy ( const Element<DataType>& element, LocalVector& lv ) const
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

    ////////////////////////////////////////////////
    //// Specific functions for Boussinesq /////////

    virtual void assemble_local_scalar_temp_mean ( const Element<DataType>& element, DataType& ls ) const
    {
        ;
    }

    virtual void assemble_local_scalar_rad_heat_volume ( const Element<DataType>& element, DataType& ls ) const
    {
        ;
    }

    virtual void assemble_local_scalar_rad_heat_surface ( const Element<DataType>& element, int facet_number, LocalVector& lv ) const
    {
        ;
    }

    void compute_tau_temp_supg ( const Element<DataType>& element, DataType& tau ) const;

    // temperature parameters
    DataType kappa_;
    DataType alpha_g_;
    DataType ref_T_;
    DataType start_T_;
    std::vector<DataType> g_;

    DataType nusselt_r_min_;
    DataType nusselt_r_max_;
    int nusselt_surface_id_;

    TempSUPGMode temp_supg_mode_;

    DataType temp_supg_gamma_;
    std::vector< DataType > temp_supg_comp_;

    int bous_scalars_;

    DataType theta_mom_buo_c_; // theta_n * g
    DataType theta_mom_buo_p_; // theta_(n-1) * g

    DataType theta_heat_diff_c_; // kappa * Laplace(theta_n)
    DataType theta_heat_diff_p_; // kappa * Laplace(theta_(n-1))
    DataType theta_heat_adv_cc_; // u_n * nabla(theta_n)
    DataType theta_heat_adv_pc_; // u_(n-1) * nabla(theta_n)
    DataType theta_heat_adv_cp_; // u_n * nabla(theta_(n-1))
    DataType theta_heat_adv_pp_; // u_(n-1) * nabla(theta_(n-1))
    DataType theta_heat_supg_c_;
    DataType theta_heat_supg_p_;

    DataType delta_mom_sou_cc_;
    DataType delta_mom_sou_nc_;
    DataType delta_mom_sou_cn_;
    DataType delta_mom_sou_nn_;
    DataType delta_mom_sou_pc_;

    DataType delta_heat_diff_c_; // kappa * Laplace(theta_n)
    DataType delta_heat_diff_n_; // kappa * Laplace(theta_(n-1))
    DataType delta_heat_adv_cc_; // u_n * nabla(theta_n)
    DataType delta_heat_adv_nc_; // u_(n-1) * nabla(theta_n)
    DataType delta_heat_adv_cn_; // u_n * nabla(theta_(n-1))
    DataType delta_heat_adv_nn_; // u_(n-1) * nabla(theta_(n-1))
    DataType delta_heat_adv_pc_; // u_(n-1) * nabla(theta_(n-1))
    DataType delta_heat_rea_cc_;
    DataType delta_heat_rea_nc_;
    DataType delta_heat_rea_cn_;
    DataType delta_heat_rea_nn_;
    DataType delta_heat_rea_pc_;
    DataType delta_heat_sou_c_;
    DataType delta_heat_sou_n_;

    using MetFlowAssembler<DIM, DataType>::num_dofs;
    using MetFlowAssembler<DIM, DataType>::dof_index;
    using MetFlowAssembler<DIM, DataType>::phi;
    using MetFlowAssembler<DIM, DataType>::grad_phi;
};

#endif
