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

#ifndef MET_FLOW_CONVDIFF_ASSEMBLER_H
#    define MET_FLOW_CONVDIFF_ASSEMBLER_H

///
/// \brief Assembler base class for scalar convection-diffusion-reaction equations
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
#    include "../../src/mesh/refined_mesh_db_view.h"
#    include "../general/goal_functional.h"
#    include "../met_flow_assembler.h"

template<int DIM, class DataType>
class MetFlowConvDiffAssembler : public virtual MetFlowAssembler<DIM, DataType>
{
  public:

    MetFlowConvDiffAssembler ( );

    ~MetFlowConvDiffAssembler ( )
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
    virtual void set_time_discretization_to_scaling ( DataType beta1, DataType beta2, DataType beta3, DataType beta4 );

    virtual void print_parameters ( );
    virtual void print_coeff ( );
    virtual DataType get_coeff ( std::string term );

    virtual void assemble_local_matrix ( const Element<DataType>& element, LocalMatrix& lm ) const;
    virtual void assemble_local_vector ( const Element<DataType>& element, LocalVector& lv ) const;
    virtual void assemble_local_scalar ( const Element<DataType>& element, DataType& ls ) const;
    virtual void assemble_local_scalar_boundary ( const Element<DataType>& element, int facet_number, LocalVector& lv ) const;

    //////////////////////////////////////////////////////
    ///// Specific functions for convection diffsuion /////

    void set_kappa ( DataType kappa )
    {
        this->kappa_ = kappa;
    }

    void set_lambda ( DataType lambda )
    {
        this->lambda_ = lambda;
    }

    void set_gamma ( DataType gamma )
    {
        this->gamma_ = gamma;
    }

    void set_vector_mode_to_std ( )
    {
        this->vector_asm_mode_ = VECTOR_STD;
    }

    void set_vector_mode_to_goal ( )
    {
        this->vector_asm_mode_ = VECTOR_GOAL;
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

    // model parameters
    DataType kappa_; // diffusion
    DataType lambda_; // reaction
    DataType gamma_; // convection

    // different modes
    VectorMode vector_asm_mode_;
    ScaMode sca_mode_;

    // time stepping coefficients
    DataType theta_diff_c_; // nu * Laplace(u_n)
    DataType theta_diff_p_; // nu * Laplace(u_(n-1))
    DataType theta_adv_cc_; // u_n * nabla(u_n)
    DataType theta_adv_pc_; // u_(n-1) * nabla(u_n)
    DataType theta_adv_cp_; // u_n * nabla(u_(n-1))
    DataType theta_adv_pp_; // u_(n-1) * nabla(u_(n-1))
    DataType theta_rea_c_; // Cor(u_n)
    DataType theta_rea_p_; // Cor(u_(n-1))
    DataType theta_sou_c_;
    DataType theta_sou_p_;

    DataType delta_diff_c_;
    DataType delta_diff_n_;
    DataType delta_adv_cc_;
    DataType delta_adv_nc_;
    DataType delta_adv_cn_;
    DataType delta_adv_nn_;
    DataType delta_adv_pc_;
    DataType delta_rea_c_;
    DataType delta_rea_n_;

    using MetFlowAssembler<DIM, DataType>::num_dofs;
    using MetFlowAssembler<DIM, DataType>::dof_index;
    using MetFlowAssembler<DIM, DataType>::phi;
    using MetFlowAssembler<DIM, DataType>::grad_phi;
};

#endif
