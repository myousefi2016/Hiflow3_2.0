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

#ifndef MET_FLOW_ESTIMATOR_ASSEMBLER_H
#    define MET_FLOW_ESTIMATOR_ASSEMBLER_H

///
/// \brief Assembler class for DWR estimator
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
#    include "adaptivity/patch_interpolation.h"
#    include "assembly/dg_assembly_assistant.h"
#    include "assembly/dg_assembly.h"
#    include "met_flow_assembler.h"

enum CellMode
{
    RESIDUAL = 0,
    TIME_JUMP = 1
};

template<int DIM, class DataType>
class MetFlowEstimatorAssembler : public virtual DGAssemblyAssistant<DIM, DataType>, public virtual MetFlowAssembler<DIM, DataType>
{
    typedef boost::function3<void,
    const mesh::Entity&, // cell
    const std::vector<Coordinate>&, // reference coordinates
    std::vector<DataType>& // values of function at the points
    > EvalFunction;
    typedef typename DGGlobalAssembler<DataType>::InterfaceSide InterfaceSide;

  public:

    MetFlowEstimatorAssembler ( );

    ~MetFlowEstimatorAssembler ( )
    {
        ;
    }

    void set_indicators ( bool primal, bool dual, bool temporal, bool spatial, bool csi );

    void set_cell_mode ( CellMode mode )
    {
        this->cell_mode_ = mode;
    }

    void set_final_time_interval ( bool flag )
    {
        this->final_time_interval_ = flag;
    }

    void set_unique_fine_space ( bool flag )
    {
        this->unique_fine_space_ = flag;
    }

    void set_est_rel_time ( int time )
    {
        this->est_rel_time_ = time;
        MetFlowAssembler<DIM, DataType>::set_rel_time ( time );
    }

    void set_time_offset ( DataType time )
    {
        this->time_offset_ = time;
    }

    void use_dwr ( bool flag )
    {
        this->use_dwr_ = flag;
    }

    void set_time_fem_order ( const std::vector<int>& order_p, const std::vector<int>& order_d )
    {
        assert ( order_p.size ( ) == this->num_var_ );
        assert ( order_d.size ( ) == this->num_var_ );
        this->time_order_P_ = order_p;
        this->time_order_D_ = order_d;
    }

    void set_fineP ( const LAD::VectorType& vector )
    {
        this->vector_fineP_ = &vector;
    }

    void set_fineP_prev ( const LAD::VectorType& vector )
    {
        this->vector_fineP_prev_ = &vector;
    }

    void set_fineP_next ( const LAD::VectorType& vector )
    {
        this->vector_fineP_next_ = &vector;
    }

    void set_fineD ( const LAD::VectorType& vector )
    {
        this->vector_fineD_ = &vector;
    }

    void set_fineD_next ( const LAD::VectorType& vector )
    {
        this->vector_fineD_next_ = &vector;
    }

    void set_fineD_prev ( const LAD::VectorType& vector )
    {
        this->vector_fineD_prev_ = &vector;
    }

    void set_fineP_space ( /*const*/ VectorSpace<DataType>& space )
    {
        this->space_fineP_ = &space;
    }

    void set_fineP_space_next ( /*const*/ VectorSpace<DataType>& space )
    {
        this->space_fineP_next_ = &space;
    }

    void set_fineP_space_prev ( /*const*/ VectorSpace<DataType>& space )
    {
        this->space_fineP_prev_ = &space;
    }

    void set_fineD_space ( /*const*/ VectorSpace<DataType>& space )
    {
        this->space_fineD_ = &space;
    }

    void set_fineD_space_next ( /*const*/ VectorSpace<DataType>& space )
    {
        this->space_fineD_next_ = &space;
    }

    void set_fineD_space_prev ( /*const*/ VectorSpace<DataType>& space )
    {
        this->space_fineD_prev_ = &space;
    }

    virtual void setup_time_quadrature ( int order );
    virtual void setup_grid_search ( );
    virtual void setup_fe_evaluators ( );

    virtual void operator() ( const Element<DataType>& element, const Quadrature<DataType>& quadrature, DataType& ls )
    {
        this->initialize_for_element_diameter ( element, quadrature );
        ls = this->h_K_;
    }

    virtual void operator() ( const Element<DataType>& element, const Quadrature<DataType>& quadrature, LocalVector& vals ) = 0;

    virtual void operator() ( const Element<DataType>& left_elem,
            const Element<DataType>& right_elem,
            const Quadrature<DataType>& left_quad,
            const Quadrature<DataType>& right_quad,
            int left_facet_number,
            int right_facet_number,
            InterfaceSide left_if_side,
            InterfaceSide right_if_side,
            LocalVector& vals ) = 0;

    virtual void clear ( );

  protected:

    virtual void initialize_for_element_diameter ( const Element<DataType>& element, const Quadrature<DataType>& quadrature );

    virtual void initialize_for_element_estimator ( const Element<DataType>& element, const Quadrature<DataType>& quadrature );
    virtual void initialize_for_interface_estimator ( const Element<DataType>& left_elem,
                                                      const Element<DataType>& right_elem,
                                                      const Quadrature<DataType>& left_quad,
                                                      const Quadrature<DataType>& right_quad,
                                                      int left_facet_number,
                                                      int right_facet_number,
                                                      InterfaceSide left_if_side,
                                                      InterfaceSide right_if_side );

    virtual void allocate_function_evaluators ( int num_var );

    virtual void assemble_local_est_cell ( const Element<DataType>& element, const Quadrature<DataType>& quadrature, LocalVector& vals );

    virtual void assemble_local_est_interface ( const Element<DataType>& left_elem,
                                                const Element<DataType>& right_elem,
                                                const Quadrature<DataType>& left_quad,
                                                const Quadrature<DataType>& right_quad,
                                                int left_facet_number,
                                                int right_facet_number,
                                                InterfaceSide left_if_side,
                                                InterfaceSide right_if_side,
                                                LocalVector& vals );

    // **************************************
    virtual void assemble_primal_cell_residual ( const Element<DataType>& element, int eq, int est_mode, DataType& r, DataType& w ) const = 0;
    virtual void assemble_dual_cell_residual ( const Element<DataType>& element, int eq, int est_mode, DataType& r, DataType& w ) const = 0;

    virtual void assemble_dual_cell_time_jump ( const Element<DataType>& element, int eq, int est_mode, DataType& r, DataType& w ) const = 0;

    virtual void assemble_primal_interface_jump ( const Element<DataType>& left_elem,
                                                  const Element<DataType>& right_elem,
                                                  const Quadrature<DataType>& left_quad,
                                                  const Quadrature<DataType>& right_quad,
                                                  int left_facet_number,
                                                  int right_facet_number,
                                                  InterfaceSide left_if_side,
                                                  InterfaceSide right_if_side,
                                                  int eq,
                                                  int est_mode,
                                                  DataType& r,
                                                  DataType& w ) = 0;

    virtual void assemble_dual_interface_jump ( const Element<DataType>& left_elem,
                                                const Element<DataType>& right_elem,
                                                const Quadrature<DataType>& left_quad,
                                                const Quadrature<DataType>& right_quad,
                                                int left_facet_number,
                                                int right_facet_number,
                                                InterfaceSide left_if_side,
                                                InterfaceSide right_if_side,
                                                int eq,
                                                int est_mode,
                                                DataType& r,
                                                DataType& w ) = 0;

    virtual void assemble_primal_interface_boundary ( const Element<DataType>& left_elem,
                                                      const Element<DataType>& right_elem,
                                                      const Quadrature<DataType>& left_quad,
                                                      const Quadrature<DataType>& right_quad,
                                                      int left_facet_number,
                                                      int right_facet_number,
                                                      InterfaceSide left_if_side,
                                                      InterfaceSide right_if_side,
                                                      int eq,
                                                      int est_mode,
                                                      DataType& r,
                                                      DataType& w ) = 0;

    virtual void assemble_dual_interface_boundary ( const Element<DataType>& left_elem,
                                                    const Element<DataType>& right_elem,
                                                    const Quadrature<DataType>& left_quad,
                                                    const Quadrature<DataType>& right_quad,
                                                    int left_facet_number,
                                                    int right_facet_number,
                                                    InterfaceSide left_if_side,
                                                    InterfaceSide right_if_side,
                                                    int eq,
                                                    int est_mode,
                                                    DataType& r,
                                                    DataType& w ) = 0;

    virtual DataType c_i ( int index ) const;
    virtual DataType w_i ( int index ) const;

    /// evaluate primal and dual weight and trial functions for temporal and spatial error indicator
    /// @param[in] c: relative time in [0,1]
    /// @param[in] q: quadrature point index
    /// @param[in] var: variable index
    virtual DataType weightD ( DataType c, int q, int var, int mode ) const;
    virtual DataType weightP ( DataType c, int q, int var, int mode ) const;

    virtual DataType weightP_tau ( DataType c, int q, int var ) const;
    virtual DataType weightD_tau ( DataType c, int q, int var ) const;
    virtual DataType weightP_h ( DataType c, int q, int var ) const;
    virtual DataType weightD_h ( DataType c, int q, int var ) const;

    virtual DataType trialP ( DataType c, int q, int var ) const;
    virtual DataType trialD ( DataType c, int q, int var ) const;
    virtual DataType dt_trialP ( DataType c, int q, int var ) const;
    virtual DataType dt_trialD ( DataType c, int q, int var ) const;
    virtual Vec<DIM, DataType> grad_trialP ( DataType c, int q, int var ) const;
    virtual Vec<DIM, DataType> grad_trialD ( DataType c, int q, int var ) const;
    virtual Mat<DIM, DIM, DataType> hess_trialP ( DataType c, int q, int var ) const;
    virtual Mat<DIM, DIM, DataType> hess_trialD ( DataType c, int q, int var ) const;

    virtual Vec<DIM, DataType> left_grad_trialP ( DataType c, int q, int var ) const;
    virtual Vec<DIM, DataType> right_grad_trialP ( DataType c, int q, int var ) const;

    virtual Vec<DIM, DataType> left_grad_trialD ( DataType c, int q, int var ) const;
    virtual Vec<DIM, DataType> right_grad_trialD ( DataType c, int q, int var ) const;

    LAD::VectorType const& vector_fineP ( ) const
    {
        return *this->vector_fineP_;
    }

    LAD::VectorType const& vector_fineP_prev ( ) const
    {
        return *this->vector_fineP_prev_;
    }

    LAD::VectorType const& vector_fineP_next ( ) const
    {
        return *this->vector_fineP_next_;
    }

    LAD::VectorType const& vector_fineD ( ) const
    {
        return *this->vector_fineD_;
    }

    LAD::VectorType const& vector_fineD_next ( ) const
    {
        return *this->vector_fineD_next_;
    }

    LAD::VectorType const& vector_fineD_prev ( ) const
    {
        return *this->vector_fineD_prev_;
    }

    DataType get_absolute_time ( DataType c ) const;

    void allocate_function_values ( int num_var );

    std::vector<DataType> quad_c_;
    std::vector<DataType> quad_w_;
    int quad_order_;

    // compute coefficients for patchwise time interpolation
    int est_rel_time_;
    DataType t_;
    //DataType dT_pc_;
    //DataType dT_cn_;
    bool final_time_interval_;
    DataType time_offset_;

    bool assemble_primal_indicator_;
    bool assemble_dual_indicator_;
    bool assemble_temporal_indicator_;
    bool assemble_spatial_indicator_;
    bool use_csi_;
    bool use_dwr_;
    CellMode cell_mode_;

    LAD::VectorType const* vector_fineP_; // t_i
    LAD::VectorType const* vector_fineP_prev_; // t_{i-1}
    LAD::VectorType const* vector_fineP_next_; // t_{i-1}
    LAD::VectorType const* vector_fineD_; // t_i
    LAD::VectorType const* vector_fineD_prev_;
    LAD::VectorType const* vector_fineD_next_; // t_{i-1}

    // for evaluating fine functions
    VectorSpace<DataType> const* space_fineP_;
    VectorSpace<DataType> const* space_fineP_next_;
    VectorSpace<DataType> const* space_fineP_prev_;
    VectorSpace<DataType> const* space_fineD_;
    VectorSpace<DataType> const* space_fineD_prev_;
    VectorSpace<DataType> const* space_fineD_next_;
    bool unique_fine_space_;

    GridGeometricSearch* search_fineP_;
    GridGeometricSearch* search_fineP_prev_;
    GridGeometricSearch* search_fineP_next_;
    GridGeometricSearch* search_fineD_;
    GridGeometricSearch* search_fineD_prev_;
    GridGeometricSearch* search_fineD_next_;

    std::vector<EvalFeFunction<LAD>*> funP_;
    std::vector<EvalFeFunction<LAD>*> funP_prev_;
    std::vector<EvalFeFunction<LAD>*> funP_next_;
    std::vector<EvalFeFunction<LAD>*> funD_;
    std::vector<EvalFeFunction<LAD>*> funD_prev_;
    std::vector<EvalFeFunction<LAD>*> funD_next_;

    int num_eq_;

    DataType h_E_;
    DataType h_K_;

    // primal higher order interpolation
    std::vector< std::vector<DataType> > fineP_;
    std::vector< std::vector<DataType> > fineP_next_;
    std::vector< std::vector<DataType> > fineP_prev_;

    // dual higher order interpolation
    std::vector< std::vector<DataType> > fineD_;
    std::vector< std::vector<DataType> > fineD_prev_;
    std::vector< std::vector<DataType> > fineD_next_;

    std::vector< FunctionValues< Vec<DIM, DataType> > > left_grad_solP_;
    std::vector< FunctionValues< Vec<DIM, DataType> > > left_grad_solP_prev_;
    std::vector< FunctionValues< Vec<DIM, DataType> > > left_grad_solP_next_;
    std::vector< FunctionValues< Vec<DIM, DataType> > > right_grad_solP_;
    std::vector< FunctionValues< Vec<DIM, DataType> > > right_grad_solP_prev_;
    std::vector< FunctionValues< Vec<DIM, DataType> > > right_grad_solP_next_;

    std::vector< FunctionValues< Vec<DIM, DataType> > > left_grad_solD_;
    std::vector< FunctionValues< Vec<DIM, DataType> > > left_grad_solD_next_;
    std::vector< FunctionValues< Vec<DIM, DataType> > > left_grad_solD_prev_;
    std::vector< FunctionValues< Vec<DIM, DataType> > > right_grad_solD_;
    std::vector< FunctionValues< Vec<DIM, DataType> > > right_grad_solD_next_;
    std::vector< FunctionValues< Vec<DIM, DataType> > > right_grad_solD_prev_;

    std::vector<int> time_order_P_;
    std::vector<int> time_order_D_;
};

#endif
