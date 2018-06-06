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

#ifndef MET_FLOW_ASSEMBLERS_H
#    define MET_FLOW_ASSEMBLERS_H

///
/// \file met_flow.h
/// \brief Assembler base classes for meteorologic application.
///
/// \author Teresa Beck, Martin Baumann, Philipp Gerstner
///

// TODO
// ---------------------------
// INNERPROD_ASM: H2
// TEHD_CYL: Vektor: OPT_AMS1,2
// INCOMP_CYL_DUAL: oseen mode
// *_CYL_DUAL: vektor: code optimierung
// *_DUAL: einschräkung an goal funktional : d_incom d_bous J != 0; d_incom d_tehd J != 0, d_tehd d_bous J != 0 nicht möglich,
// *CYL_ASSEMBLER: estimator: assemble_(primal/dual)_(mom,mass,energy,gauß)()
// Theta scheme für dual prolbem überprüfen
// Dual mode: set_time_discretization_off
// Dual Mode: Für GalerkinDC: goal_functional_final_type in assembler berücksichtigen! -> checken ob fix stimmt
// Allgemein: SEgfault bei Verwendung von hessians
// BUG: convection bei mesh change

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
#    include "../tmp_config/met_flow_vars.h"
#    include "../general/source_term.h"
#    include "../general/convection_term.h"
#    include "adaptivity/patch_interpolation.h"

enum Mode
{
    PRIMAL = 1,
    DUAL = -1,
    ESTIMATOR = 0,
    MASS = 2
};

enum ScaMode
{
    DIV_SQUARE = 1,
    DIV_MEAN = 2,
    TEMP_SQUARE = 3,
    TEMP_MEAN = 4,
    REYNOLDS_MEAN = 5,
    TEHD_FORCE = 6,
    VEL_SQUARE = 7,
    GRAD_VEL_SQUARE = 8,
    AZI_VEL_SQUARE = 9,
    RAD_VEL_SQUARE = 10,
    AXI_VEL_SQUARE = 11,
    AZI_GRAD_TEMP_SQUARE = 12,
    AZI_DEP = 13,
    RAD_DEP = 14,
    AXI_DEP = 15,
    AZI_PRESS = 16,
    RAD_PRESS = 17,
    AXI_PRESS = 18,
    AZI2RAD_CONV = 19,
    AZI_DISS = 20,
    RAD_DISS = 21,
    AXI_DISS = 22,
    AXI_GRAV = 23,
    AZI_GRADDIV = 24,
    RAD_GRADDIV = 25,
    AXI_GRADDIV = 26,
    RAD_HEAT_SURFACE = 27,
    RAD_HEAT_VOLUME = 28,
    NORM_SQUARE = 29,
    GOAL_INT = 30,
    GOAL_FIN = 31,
    INNER_PROD = 32,
    DIV_MAX = 33
};

enum GradDivMode
{
    GD_OFF = 0,
    CONST = 1,
    VMS = 2,
    SUPG = 3
};

enum TempSUPGMode
{
    TS_OFF = 0,
    STD = 1,
    DEP = 2,
    DEP_IMPLICIT = 3
};

enum SkewMode
{
    SS_OFF = 0,
    ON = 1
};

enum VectorMode
{
    VECTOR_STD = 0,
    VECTOR_VORT = 1,
    VECTOR_GOAL = 2,
    VECTOR_DEP = 3,
    VECTOR_QUANT = 4,
    VECTOR_GRADTEMP = 5,
    VECTOR_FORCE = 6
};

enum ConvectionMode
{
    STOKES = 0,
    OSEEN = 1,
    NAVIERSTOKES = 2
};

template<int DIM, class DataType>
void evaluate_general_fe_function ( const std::vector< Vec<DIM, DataType> >& quad_points,
                                    const VectorSpace<DataType>& space,
                                    const boost::function3<void, const mesh::Entity&, const std::vector<Coordinate>&, std::vector<DataType>& >& fun,
                                    SortedArray<int>& trial_cells,
                                    std::vector<DataType>& vals )
{
    int rank;
    MPI_Comm_rank ( MPI_COMM_WORLD, &rank );

    const int num_q = quad_points.size ( );
    vals.resize ( num_q, 0. );

    if ( fun == NULL )
    {
        return;
    }

    PointEvaluator<DataType> eval ( space );

    for ( int q = 0; q < num_q; ++q )
    {
        std::vector<Coordinate> point ( DIM );
        for ( int l = 0; l < DIM; ++l )
        {
            point[l] = quad_points[q][l];
        }

        DataType val = 0.;
        std::vector<int> cells;
        eval.set_trial_cells ( trial_cells.data ( ) );
        bool found = eval.evaluate_fun ( fun, point, val, cells );

        //std::cout << "EVALUATOR " << found << " " << val << " (" << point[0] << " , " << point[1] << ") CELLS: ";
        //for (int l=0; l<cells.size(); ++l)
        //{
        //    std::cout << cells[l] << " ";
        //}
        if ( !found )
        {
            std::cout << "[" << rank << "] EVALUATOR NOT FOUND " << point[0] << " , " << point[1] << std::endl;
            val = 0.;
            //found = eval.evaluate_fun_global(fun, point, val, cells, space.get_mpi_comm());
        }
        vals[q] = val;

        for ( int l = 0; l < cells.size ( ); ++l )
        {
            trial_cells.find_insert ( cells[l] );
        }
    }
}

/// \brief Abstract base class for met flow assembler implementations.
///

template<int DIM, class DataType>
class MetFlowAssembler : public virtual AssemblyAssistant<DIM, DataType>
{
    typedef boost::function3<void,
    const mesh::Entity&, // cell
    const std::vector<Coordinate>&, // reference coordinates
    std::vector<DataType>& // values of function at the points
    > EvalFunction;

    typedef boost::function4<void, int, const DataType, const Vec<DIM, DataType>&, DataType& > SourceFunction;
    // variable, time, point, return value

  public:
    MetFlowAssembler ( );

    ~MetFlowAssembler ( );

    void set_mode_to_primal ( )
    {
        mode_ = PRIMAL;
    }

    void set_mode_to_dual ( )
    {
        mode_ = DUAL;
    }

    void set_mode_to_mass ( )
    {
        mode_ = MASS;
    }

    void set_perturb_scale ( DataType scale )
    {
        this->perturb_scale_ = scale;
    }

    void set_dT_pc ( DataType dt )
    {
        this->dT_pc_ = dt;
    }

    void set_dT_cn ( DataType dt )
    {
        this->dT_cn_ = dt;
    }

    void set_time ( DataType t )
    {
        this->t_ = t;
    }

    void set_goal_functional ( GoalFunctional<DIM, DataType>* j )
    {
        goal_functional_ = j;
    }

    void set_rel_time ( int time )
    {
        this->rel_time_ = time;
    }

    void set_final_time_goal_contrib ( bool flag )
    {
        this->final_goal_contrib_ = flag;
    }

    void set_scalar_mode_to_goal_int ( )
    {
        sca_mode_ = GOAL_INT;
    }

    void set_scalar_mode_to_goal_fin ( )
    {
        sca_mode_ = GOAL_FIN;
    }

    DataType get_dT_pc ( )
    {
        return this->dT_pc_;
    }

    DataType get_dT_cn ( )
    {
        return this->dT_cn_;
    }

    void set_perturb_type ( int var, bool L2, bool H1 );

    std::vector<bool> get_L2_perturb_flags ( )
    {
        return this->L2_perturb_;
    }

    std::vector<bool> get_H1_perturb_flags ( )
    {
        return this->H1_perturb_;
    }
    void print_perturb ( );

    virtual void set_time_discretization_to_zero ( );
    virtual void set_time_discretization_to_simple ( );
    virtual void set_time_discretization_to_theta ( DataType theta );
    virtual void set_time_discretization_to_galerkin_cd ( bool modified );
    virtual void set_time_discretization_to_galerkin_dc ( bool modified );
    virtual void set_time_discretization_off ( );

    virtual void print_parameters ( );
    virtual void print_coeff ( );
    virtual DataType get_coeff ( std::string term );

    void set_source_term ( SourceTerm<DIM, DataType>* f )
    {
        this->source_term_ = f;
    };

    void set_convection_term ( ConvectionTerm<DIM, DataType>* f )
    {
        this->convection_term_ = f;
    };

    void set_newton_solution ( const LAD::VectorType& newton_sol );

    void set_solP ( const LAD::VectorType& vector )
    {
        this->vector_solP_ = &vector;
    }

    void set_solP_prev ( const LAD::VectorType& vector )
    {
        this->vector_solP_prev_ = &vector;
    }

    void set_solP_next ( const LAD::VectorType& vector )
    {
        this->vector_solP_next_ = &vector;
    }

    void set_solD ( const LAD::VectorType& vector )
    {
        this->vector_solD_ = &vector;
    }

    void set_solD_prev ( const LAD::VectorType& vector )
    {
        this->vector_solD_prev_ = &vector;
    }

    void set_solD_next ( const LAD::VectorType& vector )
    {
        this->vector_solD_next_ = &vector;
    }

    void set_perturb ( const LAD::VectorType& vector )
    {
        this->vector_perturb_ = &vector;
    }

    void set_perturb_prev ( const LAD::VectorType& vector )
    {
        this->vector_perturb_prev_ = &vector;
    }

    void set_base ( const LAD::VectorType& vector )
    {
        this->vector_base_ = &vector;
    }

    void set_conv ( const LAD::VectorType& vector )
    {
        this->vector_conv_ = &vector;
    }

    void set_conv_prev ( const LAD::VectorType& vector )
    {
        this->vector_conv_prev_ = &vector;
    }

    void set_conv_next ( const LAD::VectorType& vector )
    {
        this->vector_conv_next_ = &vector;
    }

    void set_conv_space ( /*const*/ VectorSpace<DataType>& space )
    {
        this->space_conv_ = &space;
    }

    virtual void operator() ( const Element<DataType>& element, const Quadrature<DataType>& quadrature, LocalMatrix& lm )
    {
    };

    virtual void operator() ( const Element<DataType>& element, const Quadrature<DataType>& quadrature, LocalVector& lv )
    {
    };

    virtual void operator() ( const Element<DataType>& element, const Quadrature<DataType>& quadrature, DataType& ls )
    {
    };

    virtual void operator() ( const Element<DataType>& element, int facet_number, const Quadrature<DataType>& quadrature, LocalVector& lv )
    {
    };

    virtual void operator() ( const Element<DataType>& element, const Quadrature<DataType>& quadrature, LocalVector& s_tau, LocalVector& s_h )
    {
    };

    virtual void assemble_local_matrix ( const Element<DataType>& element, LocalMatrix& lm ) const
    {
    };

    virtual void assemble_local_vector ( const Element<DataType>& element, LocalVector& lv ) const
    {
    };

    virtual void assemble_local_scalar ( const Element<DataType>& element, DataType& ls ) const
    {
    };

    virtual void assemble_local_scalar_boundary ( const Element<DataType>& element, int facet_number, LocalVector& lv ) const
    {
    };

    virtual void initialize_for_element ( const Element<DataType>& element, const Quadrature<DataType>& quadrature );
    virtual void initialize_for_element_convection ( );
    virtual void initialize_for_element_source ( );

    virtual void initialize_for_facet ( const Element<DataType>& element, const Quadrature<DataType>& quadrature, int facet_number );

    virtual void setup_grid_search ( );
    virtual void setup_fe_evaluators ( );
    virtual void setup_time_fem_evaluator ( );

    bool is_matrix_constant ( )
    {
        return is_matrix_constant_;
    }

    virtual DataType source ( DataType c, int q, int var ) const;
    virtual DataType convection ( DataType c, int q, int var ) const;
    virtual Vec<DIM, DataType> grad_convection ( DataType c, int q, int var ) const;

    virtual void clear ( );

  protected:

    void allocate_function_values ( int num_var );

    bool is_matrix_constant_;
    Mode mode_;
    ScaMode sca_mode_;

    int num_scalars_;
    int num_var_;

    DataType perturb_scale_;
    std::vector<bool> L2_perturb_;
    std::vector<bool> H1_perturb_;

    bool base_vector_set_;

    // *****************************************************
    // Time discretization stuff
    int galerkin_mode_; //(0: cont trial, disc. test; 1: discont trial, cont. test

    int rel_time_;
    DataType t_;
    DataType dT_pc_;
    DataType dT_cn_;

    DataType theta_d_dt_u_; // 0. -> stationary problem (IE-scheme!), 1. -> else
    DataType delta_d_dt_u_; // 0. -> stationary problem (IE-scheme!), 1. -> else

    DataType delta_j_c_c_;
    DataType delta_j_c_pc_;
    DataType delta_j_n_c_;
    DataType delta_j_n_nc_;
    DataType delta_j_n_n_;

    bool final_goal_contrib_;

    TimePatchInterpolation< DataType, DataType> scalar_eval_;
    TimePatchInterpolation< Vec<DIM, DataType>, DataType> vec_eval_;
    TimePatchInterpolation< Mat<DIM, DIM, DataType>, DataType> mat_eval_;

    // *****************************************************
    GoalFunctional<DIM, DataType>* goal_functional_;

    std::vector<EvalFeFunction<LAD>*> funConv_;
    std::vector<EvalFeFunction<LAD>*> funConv_prev_;
    std::vector<EvalFeFunction<LAD>*> funConv_next_;

    GridGeometricSearch* search_conv_;
    VectorSpace<DataType> const* space_conv_;

    std::vector< std::vector<EvalDerivativeFeFunction<LAD, DIM>* > > funGradConv_;
    std::vector< std::vector<EvalDerivativeFeFunction<LAD, DIM>* > > funGradConv_prev_;
    std::vector< std::vector<EvalDerivativeFeFunction<LAD, DIM>* > > funGradConv_next_;

    SourceTerm<DIM, DataType>* source_term_;
    ConvectionTerm<DIM, DataType>* convection_term_;

    // *****************************************************
    // get and set coeffcient vectors
    LAD::VectorType const* vector_solP_; // t_i
    LAD::VectorType const* vector_solP_prev_; // t_{i-1}
    LAD::VectorType const* vector_solP_next_; // t_{i+1}
    LAD::VectorType const* vector_solD_; // t_i
    LAD::VectorType const* vector_solD_next_; // t_{i+1}
    LAD::VectorType const* vector_solD_prev_; // t_{i+1}
    LAD::VectorType const* vector_perturb_; // perturbation at t_i
    LAD::VectorType const* vector_perturb_prev_; // perturbation at t_{i-1}
    LAD::VectorType const* vector_base_; // t_i
    LAD::VectorType const* vector_conv_next_;
    LAD::VectorType const* vector_conv_;
    LAD::VectorType const* vector_conv_prev_;

    LAD::VectorType const& vector_solP ( ) const
    {
        return *this->vector_solP_;
    } // t_i

    LAD::VectorType const& vector_solP_prev ( ) const
    {
        return *this->vector_solP_prev_;
    } // t_{i-1}

    LAD::VectorType const& vector_solP_next ( ) const
    {
        return *this->vector_solP_next_;
    } // t_{i+1}

    LAD::VectorType const& vector_solD ( ) const
    {
        return *this->vector_solD_;
    } // t_i

    LAD::VectorType const& vector_solD_next ( ) const
    {
        return *this->vector_solD_next_;
    } // t_{i+1}

    LAD::VectorType const& vector_solD_prev ( ) const
    {
        return *this->vector_solD_prev_;
    } // t_{i+1}

    LAD::VectorType const& vector_perturb ( ) const
    {
        return *this->vector_perturb_;
    } // t_i

    LAD::VectorType const& vector_perturb_prev ( ) const
    {
        return *this->vector_perturb_prev_;
    } // t_{i-1}

    LAD::VectorType const& vector_base ( ) const
    {
        return *this->vector_base_;
    }

    LAD::VectorType const& vector_conv ( ) const
    {
        return *this->vector_conv_;
    } // t_i

    LAD::VectorType const& vector_conv_prev ( ) const
    {
        return *this->vector_conv_prev_;
    } // t_{i-1}

    LAD::VectorType const& vector_conv_next ( ) const
    {
        return *this->vector_conv_next_;
    } // t_{i+1}

    // primal variables
    std::vector< FunctionValues<DataType> > solP_; // t_i
    std::vector< FunctionValues<DataType> > solP_prev_; // t_i-1
    std::vector< FunctionValues<DataType> > solP_next_; // t_i+1

    std::vector< FunctionValues< Vec<DIM, DataType> > > grad_solP_;
    std::vector< FunctionValues< Vec<DIM, DataType> > > grad_solP_prev_;
    std::vector< FunctionValues< Vec<DIM, DataType> > > grad_solP_next_;

    std::vector< FunctionValues< Mat<DIM, DIM, DataType> > > hess_solP_;
    std::vector< FunctionValues< Mat<DIM, DIM, DataType> > > hess_solP_prev_;
    std::vector< FunctionValues< Mat<DIM, DIM, DataType> > > hess_solP_next_;

    // dual variables
    std::vector< FunctionValues<DataType> > solD_;
    std::vector< FunctionValues<DataType> > solD_prev_;
    std::vector< FunctionValues<DataType> > solD_next_;

    std::vector< FunctionValues< Vec<DIM, DataType> > > grad_solD_;
    std::vector< FunctionValues< Vec<DIM, DataType> > > grad_solD_prev_;
    std::vector< FunctionValues< Vec<DIM, DataType> > > grad_solD_next_;

    std::vector< FunctionValues<Mat<DIM, DIM, DataType> > > hess_solD_;
    std::vector< FunctionValues<Mat<DIM, DIM, DataType> > > hess_solD_prev_;
    std::vector< FunctionValues<Mat<DIM, DIM, DataType> > > hess_solD_next_;

    // perturbation
    std::vector< FunctionValues<DataType> > perturb_;
    std::vector< FunctionValues<DataType> > perturb_prev_;

    std::vector< FunctionValues< Vec<DIM, DataType> > > grad_perturb_;
    std::vector< FunctionValues< Vec<DIM, DataType> > > grad_perturb_prev_;

    // base state
    std::vector< FunctionValues<DataType> > base_;

    // Function values
    std::vector<DataType> conv_[DIM]; // t_i
    std::vector<DataType> conv_prev_[DIM]; // t_i-1
    std::vector<DataType> conv_next_[DIM]; // t_i+1
    std::vector< Vec<DIM, DataType> > grad_conv_[DIM];
    std::vector< Vec<DIM, DataType> > grad_conv_prev_[DIM];
    std::vector< Vec<DIM, DataType> > grad_conv_next_[DIM];

    // source function
    std::vector< std::vector<DataType> > source_;
    std::vector< std::vector<DataType> > source_prev_;
    std::vector< std::vector<DataType> > source_next_;

};

#endif
