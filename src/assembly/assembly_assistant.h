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

#ifndef _ASSEMBLY_ASSISTANT_H_
#    define _ASSEMBLY_ASSISTANT_H_

#    include <algorithm>
#    include <numeric>

#    include "assembly/assembly_assistant_values.h"
#    include "assembly/function_values.h"
#    include "common/pointers.h"
#    include "common/vector_algebra.h"
#    include "fem/cell_transformation.h"
#    include "fem/fetype.h"
#    include "mesh/cell_type.h"
#    include "mesh/types.h"
#    include "quadrature/quadrature.h"
#    include "space/element.h"
#    include "linear_algebra/seq_dense_matrix.h"

///
/// \file assembly_assistant.h Convenience functionality for local assembly.
/// \author Staffan Ronnas<br>Simon Gawlok
/// \see concept_assembly

namespace hiflow
{
    using doffem::FEType;

    // TODO: if this type of assembly is desired, talk to Michael about
    // changing Quadrature, etc to use Vec & Mat classes

    /// \brief Helper class for local assembly and integration functors.
    /// \see concept_assembly

    template<int DIM, class DataType>
    class AssemblyAssistant
    {
      public:
        // Typenames for local assembly
        typedef la::SeqDenseMatrix<DataType> LocalMatrix;
        typedef std::vector<DataType> LocalVector;

        // TODO: let user decide what to compute through parameter
        AssemblyAssistant ( );

        // TODO: move these functions outside class and put documentation in here.

        /// \brief Shape function values on physical element.
        ///
        /// \param s    index of the shape function (relative to variable var if given)
        /// \param q    index of the quadrature point
        /// \param var  variable number
        /// \return     Shape function value \f$\varphi_{v,i}(\xi_q)\f$
        inline DataType phi ( int s, int q, int var = 0 ) const;

        /// \brief Shape function gradients on physical element.
        ///
        /// \param s    index of the shape function (relative to variable var if given)
        /// \param q    index of the quadrature point
        /// \param var  variable number
        /// \return     Shape function gradient \f$\nabla\varphi_{v,i}(x_q)\f$
        inline const Vec<DIM, DataType>& grad_phi ( int s, int q, int var = 0 ) const;

        /// \brief Shape function hessian on physical element.
        ///
        /// \param s    index of the shape function (relative to variable var if given)
        /// \param q    index of the quadrature point
        /// \param var  variable number
        /// \return     Shape function hessian \f$\partial_j\partial_k\varphi_{v,i}(x_q)\f$
        /// \todo   This function has not been tested yet.
        inline const Mat<DIM, DIM, DataType>& H_phi ( int s, int q, int var = 0 ) const;

        /// \brief Quadrature point on physical element.
        ///
        /// \param q    index of the quadrature point
        /// \return     Quadrature point \f$x_q\f$ on physical element
        inline const Vec<DIM, DataType>& x ( int q ) const;

        /// \brief Jacobian matrix of cell transformation.
        ///
        /// \param q    index of the quadrature point
        /// \return     Jacobian matrix \f$DF_K(\xi_q)\f$ of the cell transformation \f$F_K\f$.
        inline const Mat<DIM, DIM, DataType>& J ( int q ) const;

        /// \brief Inverse transpose of Jacobian matrix of cell transformation.
        ///
        /// \param q    index of the quadrature point
        /// \return     Inverse transpose \f$(DF_K(\xi_q))^{-T}\f$ of the jacobian matrix of the cell transformation \f$F_K\f$.
        inline const Mat<DIM, DIM, DataType>& JinvT ( int q ) const;

        /// \brief Hessian tensor of cell transformation.
        ///
        /// \param q    index of the quadrature point
        ///
        /// \return Hessian tensor \f$H_{ijk} =
        /// \partial_j\partial_k{F_i}\f$ of the cell transformation
        /// \f$F_K\f$. Index i is the position in the vector, and jk
        /// the indices for the contained matrices.
        inline const std::vector< Mat<DIM, DIM, DataType> >& H_F ( int q ) const;

        /// \brief Determinant of Jacobian matrix of cell transformation.
        ///
        /// \param q    index of the quadrature point
        /// \return     Jacobian determinant \f$det(DF_K(\xi_q))\f$ of the cell transformation \f$F_K\f$.
        inline DataType detJ ( int q ) const;

        /// \brief Jacobian matrix projected onto current facet.
        ///
        /// \pre        Assembly assistant was initialized on a facet.
        /// \param q    index of the quadrature point
        /// \return     Jf(q) = J(q) * R, where R is projection matrix corresponding to current facet.
        inline const Mat < DIM, DIM - 1, DataType > & Jf ( int q ) const;

        /// \brief Surface integration element.
        ///
        /// \pre        Assembly assistant was initialized on a facet.
        /// \param q    index of the quadrature point
        /// \return     Surface integration element ds(q) = \sqrt(det(Jf^T * Jf)) on current facet.
        inline DataType ds ( int q ) const;

        /// \brief Surface normal.
        ///
        /// \pre        Assembly assistant was initialized on a facet.
        /// \param q    index of the quadrature point
        /// \return     Surface normal n on current (reference) facet.
        inline const Vec<DIM, DataType>& n ( int q ) const;

        /// \brief The Nitsche regularization parameter used in
        /// Nitsche method.
        ///
        /// \pre        Assembly assistant was initialized on a facet.
        /// \param q    index of the quadrature point
        /// \return     Nitsche regularization parameter.
        /// \todo       Compute it in some good way.
        inline DataType nitsche_regularization ( int q ) const;

        /// \brief Quadrature point on reference element.
        ///
        /// \param q    index of the quadrature point
        /// \return     Quadrature point \f$\xi_q\f$ on reference element
        inline const Vec<DIM, DataType>& q_point ( int q ) const;

        /// \brief Number of quadrature points.
        ///
        /// \return Number of quadrature points for current quadrature rule.
        inline int num_quadrature_points ( ) const;

        /// \brief Quadrature weights.
        ///
        /// \param q    index of the quadrature point
        /// \return     Quadrature weight \f$w_q\f$ for the
        inline DataType w ( int q ) const;

        /// \brief Number of local dofs for a variable.
        ///
        /// \param var  variable number
        /// \return Number of local dofs for variable var.
        inline int num_dofs ( int var ) const;

        /// \brief Number of local dofs.
        ///
        /// \return Total number of local dofs for all variables.
        inline int num_dofs_total ( ) const;

        /// \brief Number of variables.
        /// \return Number of variables.
        inline int num_vars ( ) const;

        /// \brief Local dof index for dof s of variable var.
        ///
        /// \param s     dof index relative to variable var ( in range [0, num_dofs(var)) )
        /// \param var   variable number
        /// \return Local dof index \f$l = l(var, s)\f$ for dof \f$s\f$ for variable var.
        inline int dof_index ( int s, int var ) const;

        /// \brief Element size parameter h
        /// \return Element size parameter h
        inline DataType h ( ) const;

        /// \brief Evaluate a finite element function on the current element.
        ///
        /// \param[in] coefficients  global vector of dof values for the function.
        /// \param[in] var           variable to compute the function for.
        /// \param[out] function_values  function values evaluated at all quadrature points.
        template<class VectorType>
        void evaluate_fe_function ( const VectorType& coefficients, int var,
                                    FunctionValues<DataType>& function_values ) const;

        /// \brief Evaluate the gradient of a finite element function on the current element.
        ///
        /// \param[in] coefficients  global vector of dof values for the function.
        /// \param[in] var           variable to compute the gradient for.
        /// \param[out] function_gradients  function gradients at all quadrature points.
        template<class VectorType>
        void evaluate_fe_function_gradients ( const VectorType& coefficients, int var,
                                              FunctionValues< Vec<DIM, DataType> >& function_gradients ) const;

        /// \brief Evaluate the hessian of a finite element function on the current element.
        ///
        /// \param[in] coefficients  global vector of dof values for the function.
        /// \param[in] var           variable to compute the gradient for.
        /// \param[out] function_hessians  function gradients at all quadrature points.
        template<class VectorType>
        void evaluate_fe_function_hessians ( const VectorType& coefficients, int var,
                                             FunctionValues< Mat<DIM, DIM, DataType> >& function_hessians ) const;

        /// \brief Initialize the assistant for an element and a quadrature.
        ///
        /// \details Recomputes all necessary values on the given
        /// element for the given quadrature formula
        ///
        /// \param[in] element     element for which to initialize the AssemblyAssistant
        /// \param[in] quadrature  quadrature rule
        void initialize_for_element ( const Element<DataType>& element,
                                      const Quadrature<DataType>& element_quadrature );

        /// \brief Initialize the assistant for a facet of the element and a quadrature.
        ///
        /// \details Recomputes all necessary values on the given
        /// element for the given quadrature formula and facet.
        ///
        /// \param[in] element          element for which to initialize the AssemblyAssistant
        /// \param[in] facet_quadrature quadrature rule (should have points only on the given element facet)
        /// \param[in] facet_number     local number of facet in element
        void initialize_for_facet ( const Element<DataType>& element,
                                    const Quadrature<DataType>& facet_quadrature,
                                    int facet_number );

      private:
        // initialization sub-functions
        bool initialize_quadrature ( const Quadrature<DataType>& new_quadrature, bool force_update );
        bool initialize_fe_types ( const Element<DataType>& element, bool force_update );
        bool initialize_cell_values ( const Element<DataType>& element, bool force_update );
        bool initialize_facet_values ( const Element<DataType>& element, int facet, bool force_update );

        int find_fe_type ( const FEType<DataType>* fe_type ) const;
        void compute_quadrature_values ( );
        void compute_fe_values ( );
        void compute_transformed_cell_values ( );
        void compute_transformed_facet_values ( );
        void compute_facet_projection_matrix ( );
        void compute_facet_normal ( );

        template<class VectorType>
        std::vector<DataType> extract_dof_values ( const VectorType& coefficients, int var ) const;

        mesh::EntityNumber cell_index_; ///< Current cell index
        int facet_number_; ///< Current facet number (-1 on cell)
        Mat < DIM, DIM - 1, DataType > facet_proj_; ///< Current facet projection matrix
        Vec<DIM, DataType> n_; ///< Current facet normal
        FunctionValues< Vec<DIM, DataType> > mapped_n_; ///< Normals on physical element
        DataType nitsche_regularization_; ///< Regularization parameter for Nitsche method
        Quadrature<DataType> quadrature_; ///< Current quadrature
        std::vector<const FEType<DataType>* > fe_types_; ///< Current FEType:s
        std::vector<int> fe_offsets_for_var_; ///< Mapping var -> FEType
        std::vector<int> dof_offset_for_var_; ///< Mapping var -> dof offset
        std::vector<int> num_dofs_for_var_; ///< Number of dofs / var
        std::vector<int> global_dof_indices_; ///< Global dofs for current element

        const CellTransformation<DataType>* cell_transform_; ///< Cell transform for current element
        mesh::CellType::Tag cell_type_; ///< Cell type of current element
        DataType h_; ///< Element size parameter h

        std::vector< Vec<DIM, DataType> > q_points_; ///< Quadrature points on reference element
        std::vector<DataType> weights_; ///< Quadrature weights

        FunctionValues< std::vector< DataType> > phi_; ///< Shape function values

        ///< Shape function gradients on reference element
        FunctionValues< std::vector< Vec<DIM, DataType> > > grad_phi_hat_;

        ///< Shape function hessians on reference element
        FunctionValues< std::vector< Mat<DIM, DIM, DataType> > > H_phi_hat_;

        ///< Shape function gradients on physical element
        FunctionValues< std::vector< Vec<DIM, DataType> > > grad_phi_;
        FunctionValues< Vec<DIM, DataType> > x_; ///< Quadrature points on physical element
        FunctionValues< Mat<DIM, DIM, DataType> > J_; ///< Jacobian matrix at quadrature points
        FunctionValues<DataType> detJ_; ///< Determinants of jacobian matrices
        FunctionValues< Mat<DIM, DIM, DataType> > JinvT_; ///< Inverse-transpose of jacobian matrices
        ///< Hessian of cell transformation at quadrature points
        FunctionValues< std::vector< Mat<DIM, DIM, DataType> > > H_transform_;

        ///< H-mapped shape function gradients on physical element
        FunctionValues< std::vector< Mat<DIM, DIM, DataType> > > H_mapped_grad_;

        ///< Shape function hessians on physical element
        FunctionValues< std::vector< Mat<DIM, DIM, DataType> > > H_phi_;

        // Facet-specific values
        FunctionValues < Mat < DIM, DIM - 1, DataType > > Jf_; ///< Jacobian matrix projected onto facet
        FunctionValues<DataType> ds_; ///< Surface element on facet
    };

    template<int DIM, class DataType>
    AssemblyAssistant<DIM, DataType>::AssemblyAssistant ( )
    : cell_transform_ ( 0 )
    {
        assert ( DIM > 0 );
        assert ( DIM <= 3 );
        facet_number_ = -1;
    }

    template<int DIM, class DataType>
    inline DataType AssemblyAssistant<DIM, DataType>::phi ( int s, int q, int var ) const
    {
        return phi_[q][fe_offsets_for_var_[var] + s];
    }

    template<int DIM, class DataType>
    inline const Vec<DIM, DataType>& AssemblyAssistant<DIM, DataType>::grad_phi ( int s, int q, int var ) const
    {
        return grad_phi_[q][fe_offsets_for_var_[var] + s];
    }

    template<int DIM, class DataType>
    inline const Mat<DIM, DIM, DataType>& AssemblyAssistant<DIM, DataType>::H_phi ( int s, int q, int var ) const
    {
        return H_phi_[q][fe_offsets_for_var_[var] + s];
    }

    template<int DIM, class DataType>
    inline const Vec<DIM, DataType>& AssemblyAssistant<DIM, DataType>::x ( int q ) const
    {
        return x_[q];
    }

    template<int DIM, class DataType>
    inline const Mat<DIM, DIM, DataType>& AssemblyAssistant<DIM, DataType>::J ( int q ) const
    {
        return J_[q];
    }

    template<int DIM, class DataType>
    inline const Mat<DIM, DIM, DataType>& AssemblyAssistant<DIM, DataType>::JinvT ( int q ) const
    {
        return JinvT_[q];
    }

    template<int DIM, class DataType>
    inline const std::vector< Mat<DIM, DIM, DataType> >& AssemblyAssistant<DIM, DataType>::H_F ( int q ) const
    {
        return H_transform_[q];
    }

    template<int DIM, class DataType>
    inline DataType AssemblyAssistant<DIM, DataType>::detJ ( int q ) const
    {
        return detJ_[q];
    }

    template<int DIM, class DataType>
    inline const Mat < DIM, DIM - 1, DataType > & AssemblyAssistant<DIM, DataType>::Jf ( int q ) const
    {
        assert ( facet_number_ > -1 );
        return Jf_[q];
    }

    template<int DIM, class DataType>
    inline DataType AssemblyAssistant<DIM, DataType>::ds ( int q ) const
    {
        assert ( facet_number_ > -1 );
        return ds_[q];
    }

    template<int DIM, class DataType>
    inline const Vec<DIM, DataType>& AssemblyAssistant<DIM, DataType>::n ( int q ) const
    {
        assert ( facet_number_ > -1 );
        return mapped_n_[q];
    }

    template<int DIM, class DataType>
    inline DataType AssemblyAssistant<DIM, DataType>::nitsche_regularization ( int q ) const
    {
        assert ( facet_number_ > -1 );
        // TODO: use q and specify this parameter somehow
        return nitsche_regularization_;
    }

    template<int DIM, class DataType>
    inline const Vec<DIM, DataType>& AssemblyAssistant<DIM, DataType>::q_point ( int q ) const
    {
        return q_points_[q];
    }

    template<int DIM, class DataType>
    inline int AssemblyAssistant<DIM, DataType>::num_quadrature_points ( ) const
    {
        return q_points_.size ( );
    }

    template<int DIM, class DataType>
    inline DataType AssemblyAssistant<DIM, DataType>::w ( int q ) const
    {
        return weights_[q];
    }

    template<int DIM, class DataType>
    inline int AssemblyAssistant<DIM, DataType>::num_dofs ( int var ) const
    {
        return num_dofs_for_var_[var];
    }

    template<int DIM, class DataType>
    inline int AssemblyAssistant<DIM, DataType>::num_dofs_total ( ) const
    {
        return std::accumulate ( num_dofs_for_var_.begin ( ), num_dofs_for_var_.end ( ), 0 );
    }

    template<int DIM, class DataType>
    inline int AssemblyAssistant<DIM, DataType>::num_vars ( ) const
    {
        return num_dofs_for_var_.size ( );
    }

    template<int DIM, class DataType>
    inline int AssemblyAssistant<DIM, DataType>::dof_index ( int i, int var ) const
    {
        return dof_offset_for_var_[var] + i;
    }

    template<int DIM, class DataType>
    inline DataType AssemblyAssistant<DIM, DataType>::h ( ) const
    {
        return h_;
    }

    template<int DIM, class DataType>
    template<class VectorType>
    void AssemblyAssistant<DIM, DataType>::evaluate_fe_function ( const VectorType& coefficients, int var,
                                                                  FunctionValues<DataType>& function_values ) const
    {
        // extract dof-values corresponding to local variable
        const std::vector<DataType> local_coefficients = extract_dof_values ( coefficients, var );

        // compute function values
        function_values.compute ( phi_,
                                  EvalFiniteElementFunction<DataType>( fe_offsets_for_var_[var],
                                  local_coefficients ) );
    }

    template<int DIM, class DataType>
    template<class VectorType>
    void AssemblyAssistant<DIM, DataType>::evaluate_fe_function_gradients ( const VectorType& coefficients, int var,
                                                                            FunctionValues< Vec<DIM, DataType> >& function_gradients ) const
    {
        // extract dof-values corresponding to local variable
        const std::vector<DataType> local_coefficients = extract_dof_values ( coefficients, var );

        // compute function values
        function_gradients.compute ( grad_phi_,
                                     EvalFiniteElementFunctionGradient<DIM, DataType>( fe_offsets_for_var_[var],
                                     local_coefficients ) );
    }

    template<int DIM, class DataType>
    template<class VectorType>
    void AssemblyAssistant<DIM, DataType>::evaluate_fe_function_hessians ( const VectorType& coefficients, int var,
                                                                           FunctionValues< Mat<DIM, DIM, DataType> >& function_hessians ) const
    {
        // extract dof-values corresponding to local variable
        const std::vector<DataType> local_coefficients = extract_dof_values ( coefficients, var );

        // compute function values
        function_hessians.compute ( H_phi_,
                                    EvalFiniteElementFunctionHessian<DIM, DataType>( fe_offsets_for_var_[var],
                                    local_coefficients ) );
    }

    template<int DIM, class DataType>
    void AssemblyAssistant<DIM, DataType>::initialize_for_element ( const Element<DataType>& element,
                                                                    const Quadrature<DataType>& element_quadrature )
    {

        facet_number_ = -1;

        // initialize quadrature
        const bool changed_quadrature = initialize_quadrature ( element_quadrature, false );

        // initialize fe types
        const bool changed_fe_types = initialize_fe_types ( element, changed_quadrature );

        // compute transform values for new cell
        initialize_cell_values ( element, changed_quadrature || changed_fe_types );

        // update global dof indices
        this->global_dof_indices_.clear ( );
        element.get_dof_indices ( this->global_dof_indices_ );

    }

    template<int DIM, class DataType>
    void AssemblyAssistant<DIM, DataType>::initialize_for_facet ( const Element<DataType>& element,
                                                                  const Quadrature<DataType>& facet_quadrature,
                                                                  int facet_number )
    {

        // initialize quadrature
        bool force_update = true; //(facet_number_ != facet_number);
        facet_number_ = facet_number;

        const bool changed_quadrature = initialize_quadrature ( facet_quadrature, force_update );

        // initialize fe types
        const bool changed_fe_types = initialize_fe_types ( element, changed_quadrature );

        // compute transform values for new facet
        initialize_facet_values ( element, facet_number, changed_quadrature || changed_fe_types );

        // update global dof indices
        this->global_dof_indices_.clear ( );
        element.get_dof_indices ( this->global_dof_indices_ );
    }

    template<int DIM, class DataType>
    bool AssemblyAssistant<DIM, DataType>::initialize_quadrature ( const Quadrature<DataType>& new_quadrature,
                                                                   bool force_update )
    {
        // check if quadrature changed
        bool need_quadrature_update = force_update || quadrature_.size ( ) == 0;
        if ( !need_quadrature_update )
        {
            if ( new_quadrature.size ( ) != quadrature_.size ( )
                 || new_quadrature.name ( ) != quadrature_.name ( ) )
            {
                need_quadrature_update = true;
            }
        }

        if ( need_quadrature_update )
        {
            // copy quadrature
            quadrature_ = new_quadrature;
            compute_quadrature_values ( );
        }
        return need_quadrature_update;
    }

    template<int DIM, class DataType>
    bool AssemblyAssistant<DIM, DataType>::initialize_fe_types ( const Element<DataType>& element, bool force_update )
    {
        const size_t num_vars = element.get_num_variables ( );

        // fe_types_.empty() is true the first time function is called
        bool need_fe_update = force_update || fe_types_.empty ( );

        // check if element type changed
        if ( !need_fe_update )
        {
            for ( size_t var = 0; var != num_vars; ++var )
            {
                const FEType<DataType>* element_fe_type = element.get_fe_type ( var );

                // check that we have Lagrange element (so that we can get degree)
                assert ( element_fe_type != 0 );

                // check if element_fe_type already exists
                const int pos = find_fe_type ( element_fe_type );
                if ( pos < 0 )
                {
                    need_fe_update = true;
                    break;
                }
            }
        }

        if ( need_fe_update )
        {
            fe_types_.clear ( );
            fe_offsets_for_var_.clear ( );
            num_dofs_for_var_.clear ( );
            dof_offset_for_var_.clear ( );

            fe_types_.reserve ( num_vars );
            fe_offsets_for_var_.reserve ( num_vars );
            dof_offset_for_var_.reserve ( num_vars );

            int offset = 0;
            int num_dofs = 0;

            std::vector<int> pos2offset;
            pos2offset.reserve ( num_vars );

            // rebuild fe_types_ and var_offsets_ vectors
            for ( size_t var = 0; var != num_vars; ++var )
            {
                const FEType<DataType>* element_fe_type = element.get_fe_type ( var );
                assert ( element_fe_type != 0 );
                const int pos = find_fe_type ( element_fe_type );
                if ( pos < 0 )
                { // this fe_type has not been added before
                    fe_types_.push_back ( element_fe_type );
                    fe_offsets_for_var_.push_back ( offset );
                    pos2offset.push_back ( offset );
                    offset += element.get_num_dofs ( var );
                }
                else
                {
                    fe_offsets_for_var_.push_back ( pos2offset[pos] );
                }
                num_dofs_for_var_.push_back ( element.get_num_dofs ( var ) );
                dof_offset_for_var_.push_back ( num_dofs );
                num_dofs += num_dofs_for_var_.back ( );
            }
            compute_fe_values ( );
        }
        return need_fe_update;
    }

    template<int DIM, class DataType>
    bool AssemblyAssistant<DIM, DataType>::initialize_cell_values ( const Element<DataType>& element,
                                                                    bool force_update )
    {
        mesh::EntityNumber new_cell_index = element.get_cell_index ( );

        if ( force_update || ( cell_index_ != new_cell_index ) )
        {
            cell_index_ = new_cell_index;
            cell_type_ = element.get_cell ( ).cell_type ( ).tag ( );
            cell_transform_ = element.get_cell_transformation ( );

            compute_transformed_cell_values ( );
            return true;
        }
        return false;
    }

    template<int DIM, class DataType>
    bool AssemblyAssistant<DIM, DataType>::initialize_facet_values ( const Element<DataType>& element,
                                                                     int facet,
                                                                     bool force_update )
    {

        if ( DIM == 1 )
        {
            throw "Facet integration not possible for DIM = 1\n";
        }

        mesh::EntityNumber new_cell_index = element.get_cell_index ( );

        if ( force_update
             || facet_number_ != facet
             || ( cell_index_ != new_cell_index ) )
        {

            //compute map: local facet dof to local cell dof

            cell_index_ = new_cell_index;
            cell_type_ = element.get_cell ( ).cell_type ( ).tag ( );
            cell_transform_ = element.get_cell_transformation ( );

            compute_transformed_cell_values ( );

            // facet - specific initialization
            facet_number_ = facet;
            // TODO: Compute this parameter in some way
            nitsche_regularization_ = 1.0;
            compute_facet_projection_matrix ( );
            compute_transformed_facet_values ( );
            compute_facet_normal ( );

            return true;
        }
        return false;
    }

    template<int DIM, class DataType>
    int AssemblyAssistant<DIM, DataType>::find_fe_type ( const FEType<DataType>* fe_type ) const
    {
        assert ( fe_type != 0 );
        int pos = 0;
        bool fe_type_found = false;
        typedef typename std::vector<const FEType<DataType>* >::const_iterator Iterator;
        for ( Iterator it = fe_types_.begin ( ), e_it = fe_types_.end ( ); it != e_it; ++it )
        {
            if ( ( *it )->get_global_id ( ) == fe_type->get_global_id ( ) )
            {
                fe_type_found = true;
                break;
            }
            ++pos;
        }

        if ( fe_type_found )
        {
            return pos;
        }
        else
        {
            return -1;
        }
    }

    template<int DIM, class DataType>
    void AssemblyAssistant<DIM, DataType>::compute_quadrature_values ( )
    {
        size_t num_q = quadrature_.size ( );
        q_points_.resize ( num_q );
        weights_.resize ( num_q );
        for ( size_t q = 0; q != num_q; ++q )
        {
            for ( size_t c = 0; c != DIM; ++c )
            {
                switch ( c )
                {
                    case 0: q_points_[q][0] = quadrature_.x ( q );
                        break;
                    case 1: q_points_[q][1] = quadrature_.y ( q );
                        break;
                    case 2: q_points_[q][2] = quadrature_.z ( q );
                        break;
                    default: assert ( false );
                        break;
                }
            }
            weights_[q] = quadrature_.w ( q );
        }
    }

    template<int DIM, class DataType>
    void AssemblyAssistant<DIM, DataType>::compute_fe_values ( )
    {
        assert ( q_points_.size ( ) == quadrature_.size ( ) );
        assert ( !fe_types_.empty ( ) );

        // compute shape function values
        phi_.compute ( q_points_, EvalShapeFunctions<DIM, DataType>( fe_types_ ) );

        // compute shape function gradients
        grad_phi_hat_.compute ( q_points_, EvalShapeFunctionGradients<DIM, DataType>( fe_types_ ) );

        // compute shape function hessians
        H_phi_hat_.compute ( q_points_, EvalShapeFunctionHessians<DIM, DataType>( fe_types_ ) );
    }

    template<int DIM, class DataType>
    void AssemblyAssistant<DIM, DataType>::compute_transformed_cell_values ( )
    {
        assert ( q_points_.size ( ) == quadrature_.size ( ) );
        assert ( cell_transform_ != 0 );

        // compute quadrature points on physical cell
        x_.compute ( q_points_, EvalPhysicalPoint<DIM, DataType>( *cell_transform_ ) );

        // compute jacobians on physical cell
        J_.compute ( q_points_, EvalPhysicalJacobian<DIM, DataType>( *cell_transform_ ) );

        // compute hessians on physical cell
        H_transform_.compute ( q_points_, EvalCellTransformationHessian<DIM, DataType>( *cell_transform_ ) );
        // compute determinant of jacobians
        detJ_.compute ( J_, EvalDeterminant<DIM, DataType>( ) );

        // compute inverse-transpose of jacobians
        JinvT_.compute ( J_, EvalInvTranspose<DIM, DataType>( ) );

        // compute JinvT * grad_phi_hat
        grad_phi_.compute ( JinvT_, EvalMappedShapeFunctionGradients<DIM, DataType>( grad_phi_hat_ ) );

        H_mapped_grad_.compute ( H_transform_, EvalHMappedGradients<DIM, DataType>( grad_phi_ ) );
        H_phi_.compute ( JinvT_, EvalMappedShapeFunctionHessians<DIM, DataType>( H_phi_hat_, H_mapped_grad_ ) );

        // compute mesh parameter h
        h_ = std::pow ( static_cast < DataType > ( std::abs ( detJ ( 0 ) ) ), static_cast < DataType > ( 1. / static_cast < DataType > ( DIM ) ) );
    }

    template<int DIM, class DataType>
    void AssemblyAssistant<DIM, DataType>::compute_transformed_facet_values ( )
    {
        assert ( q_points_.size ( ) == quadrature_.size ( ) );
        assert ( cell_transform_ != 0 );

        Jf_.compute ( J_, EvalRightMatrixMult < DIM, DIM, DIM - 1, DataType > ( facet_proj_ ) );
        ds_.compute ( Jf_, EvalSurfaceElement<DIM, DataType>( ) );

        // compute mesh parameter h
        h_ = std::pow ( static_cast < DataType > ( std::abs ( detJ ( 0 ) ) ), static_cast < DataType > ( 1. / static_cast < DataType > ( DIM ) ) );
    }

    template<int DIM, class DataType>
    void AssemblyAssistant<DIM, DataType>::compute_facet_projection_matrix ( )
    {
        using mesh::CellType;

        switch ( cell_type_ )
        {
            case CellType::TRIANGLE:
                switch ( facet_number_ )
                {
                    case 0: //bottom edge
                        facet_proj_ ( 0, 0 ) = 1.;
                        facet_proj_ ( 1, 0 ) = 0.;
                        break;
                    case 1: // diagonal
                        facet_proj_ ( 0, 0 ) = 1.;
                        facet_proj_ ( 1, 0 ) = -1.;
                        break;
                    case 2: // left edge
                        facet_proj_ ( 0, 0 ) = 0.;
                        facet_proj_ ( 1, 0 ) = 1.;
                        break;
                    default:
                        assert ( 0 );
                }
                break;
            case CellType::QUADRILATERAL:
                switch ( facet_number_ )
                {
                    case 0: // bottom edge
                    case 2: // top edge
                        facet_proj_ ( 0, 0 ) = 1.;
                        facet_proj_ ( 1, 0 ) = 0.;
                        break;
                    case 1: // right edge
                    case 3: // left edge
                        facet_proj_ ( 0, 0 ) = 0.;
                        facet_proj_ ( 1, 0 ) = 1.;
                        break;
                    default:
                        assert ( 0 );
                }
                break;
            case CellType::HEXAHEDRON:
                switch ( facet_number_ )
                {
                    case 0: // bottom
                    case 5: // top
                        facet_proj_ ( 0, 0 ) = 1.;
                        facet_proj_ ( 1, 0 ) = 0.;
                        facet_proj_ ( 2, 0 ) = 0.;

                        facet_proj_ ( 0, 1 ) = 0.;
                        facet_proj_ ( 1, 1 ) = 1.;
                        facet_proj_ ( 2, 1 ) = 0.;
                        break;

                    case 1: // front 1
                    case 4: // back  4
                        facet_proj_ ( 0, 0 ) = 1.;
                        facet_proj_ ( 1, 0 ) = 0.;
                        facet_proj_ ( 2, 0 ) = 0.;

                        facet_proj_ ( 0, 1 ) = 0.;
                        facet_proj_ ( 1, 1 ) = 0.;
                        facet_proj_ ( 2, 1 ) = 1.;
                        break;

                    case 2: // left  2
                    case 3: // right 3
                        facet_proj_ ( 0, 0 ) = 0.;
                        facet_proj_ ( 1, 0 ) = 1.;
                        facet_proj_ ( 2, 0 ) = 0.;

                        facet_proj_ ( 0, 1 ) = 0.;
                        facet_proj_ ( 1, 1 ) = 0.;
                        facet_proj_ ( 2, 1 ) = 1.;
                        break;

                    default:
                        assert ( 0 );
                }
                break;
            case CellType::TETRAHEDRON:
                switch ( facet_number_ )
                {
                    case 0: // bottom
                        facet_proj_ ( 0, 0 ) = 1.;
                        facet_proj_ ( 1, 0 ) = 0.;
                        facet_proj_ ( 2, 0 ) = 0.;

                        facet_proj_ ( 0, 1 ) = 0.;
                        facet_proj_ ( 1, 1 ) = 1.;
                        facet_proj_ ( 2, 1 ) = 0.;
                        break;

                    case 1: // front
                        facet_proj_ ( 0, 0 ) = 1.;
                        facet_proj_ ( 1, 0 ) = 0.;
                        facet_proj_ ( 2, 0 ) = 0.;

                        facet_proj_ ( 0, 1 ) = 0.;
                        facet_proj_ ( 1, 1 ) = 0.;
                        facet_proj_ ( 2, 1 ) = 1.;
                        break;

                    case 2: // left
                        facet_proj_ ( 0, 0 ) = 0.;
                        facet_proj_ ( 1, 0 ) = 1.;
                        facet_proj_ ( 2, 0 ) = 0.;

                        facet_proj_ ( 0, 1 ) = 0.;
                        facet_proj_ ( 1, 1 ) = 0.;
                        facet_proj_ ( 2, 1 ) = 1.;
                        break;

                    case 3: // back
                        /*facet_proj_(0, 0) = 1.;
                        facet_proj_(1, 0) = -1.;
                        facet_proj_(2, 0) = 0.;

                        facet_proj_(0, 1) = 0.;
                        facet_proj_(1, 1) = 1.;
                        facet_proj_(2, 1) = -1.;*/
                        facet_proj_ ( 0, 0 ) = 1.;
                        facet_proj_ ( 1, 0 ) = 0.;
                        facet_proj_ ( 2, 0 ) = -1.;

                        facet_proj_ ( 0, 1 ) = 0.;
                        facet_proj_ ( 1, 1 ) = 1.;
                        facet_proj_ ( 2, 1 ) = -1.;
                        break;

                    default:
                        assert ( 0 );
                }
                break;
            case CellType::PYRAMID:
                switch ( facet_number_ )
                {
                    case 0: // bottom
                        facet_proj_ ( 0, 0 ) = 1.;
                        facet_proj_ ( 1, 0 ) = 0.;
                        facet_proj_ ( 2, 0 ) = 0.;

                        facet_proj_ ( 1, 1 ) = 0.;
                        facet_proj_ ( 1, 1 ) = 1.;
                        facet_proj_ ( 2, 1 ) = 0.;
                        break;

                    case 1: //front
                        facet_proj_ ( 0, 0 ) = 1.;
                        facet_proj_ ( 1, 0 ) = 0.;
                        facet_proj_ ( 2, 0 ) = 0.;

                        facet_proj_ ( 0, 1 ) = 0.5;
                        facet_proj_ ( 1, 1 ) = 0.5;
                        facet_proj_ ( 2, 1 ) = 1.;
                        break;

                    case 2: //right
                        facet_proj_ ( 0, 0 ) = 0.;
                        facet_proj_ ( 1, 0 ) = 1.;
                        facet_proj_ ( 2, 0 ) = 0.;

                        facet_proj_ ( 0, 1 ) = -0.5;
                        facet_proj_ ( 1, 1 ) = 0.5;
                        facet_proj_ ( 2, 1 ) = 1.;
                        break;

                    case 3: //back
                        facet_proj_ ( 0, 0 ) = -1.;
                        facet_proj_ ( 1, 0 ) = 0.;
                        facet_proj_ ( 2, 0 ) = 0.;

                        facet_proj_ ( 0, 1 ) = -0.5;
                        facet_proj_ ( 1, 1 ) = -0.;
                        facet_proj_ ( 2, 1 ) = 1.;
                        break;

                    case 4: //left
                        facet_proj_ ( 0, 0 ) = 0.;
                        facet_proj_ ( 1, 0 ) = -1.;
                        facet_proj_ ( 2, 0 ) = 0.;

                        facet_proj_ ( 0, 1 ) = 0.5;
                        facet_proj_ ( 1, 1 ) = -0.5;
                        facet_proj_ ( 2, 1 ) = 1.;
                        break;

                    default:
                        assert ( 0 );

                }
                break;

            default:
                // not yet supported
                std::cerr << "Failed: Cell type " << cell_type_ << " not yet implemented!\n";
                assert ( false );
        }
    }

    template<int DIM, class DataType>
    void AssemblyAssistant<DIM, DataType>::compute_facet_normal ( )
    {
        using mesh::CellType;

        switch ( cell_type_ )
        {
            case CellType::TRIANGLE:
                switch ( facet_number_ )
                {
                    case 0: // bottom edge
                        for ( size_t d = 0; d != DIM; ++d )
                        {
                            switch ( d )
                            {
                                case 0: n_[d] = 0.0;
                                    break;
                                case 1: n_[d] = -1.0;
                                    break;
                                case 2: n_[d] = 0.0;
                                    break;
                                default: assert ( false );
                                    break;
                            }
                        }
                        break;
                    case 1: // diagonal
                        for ( size_t d = 0; d != DIM; ++d )
                        {
                            switch ( d )
                            {
                                case 0: n_[d] = 1. / std::sqrt ( 2. );
                                    break;
                                case 1: n_[d] = 1. / std::sqrt ( 2. );
                                    break;
                                case 2: n_[d] = .0;
                                    break;
                                default: assert ( false );
                                    break;
                            }
                        }
                        break;
                    case 2: // left edge
                        for ( size_t d = 0; d != DIM; ++d )
                        {
                            switch ( d )
                            {
                                case 0: n_[d] = -1.0;
                                    break;
                                case 1: n_[d] = 0.0;
                                    break;
                                case 2: n_[d] = 0.0;
                                    break;
                                default: assert ( false );
                                    break;
                            }
                        }
                        break;
                    default:
                        assert ( 0 );
                }
                break;
            case CellType::QUADRILATERAL:
                switch ( facet_number_ )
                {
                    case 0: // bottom edge
                        for ( size_t d = 0; d != DIM; ++d )
                        {
                            switch ( d )
                            {
                                case 0: n_[d] = 0.0;
                                    break;
                                case 1: n_[d] = -1.0;
                                    break;
                                case 2: n_[d] = 0.0;
                                    break;
                                default: assert ( false );
                                    break;
                            }
                        }
                        break;
                    case 1: // right edge
                        for ( size_t d = 0; d != DIM; ++d )
                        {
                            switch ( d )
                            {
                                case 0: n_[d] = 1.0;
                                    break;
                                case 1: n_[d] = 0.0;
                                    break;
                                case 2: n_[d] = 0.0;
                                    break;
                                default: assert ( false );
                                    break;
                            }
                        }
                        break;
                    case 2: // top edge
                        for ( size_t d = 0; d != DIM; ++d )
                        {
                            switch ( d )
                            {
                                case 0: n_[d] = 0.0;
                                    break;
                                case 1: n_[d] = 1.0;
                                    break;
                                case 2: n_[d] = 0.0;
                                    break;
                                default: assert ( false );
                                    break;
                            }
                        }
                        break;
                    case 3: // left edge
                        for ( size_t d = 0; d != DIM; ++d )
                        {
                            switch ( d )
                            {
                                case 0: n_[d] = -1.0;
                                    break;
                                case 1: n_[d] = 0.0;
                                    break;
                                case 2: n_[d] = 0.0;
                                    break;
                                default: assert ( false );
                                    break;
                            }
                        }
                        break;
                    default:
                        assert ( 0 );
                }
                break;

            case CellType::HEXAHEDRON:
                switch ( facet_number_ )
                {
                    case 0: // bottom
                        for ( size_t d = 0; d != DIM; ++d )
                        {
                            switch ( d )
                            {
                                case 0: n_[d] = 0.0;
                                    break;
                                case 1: n_[d] = 0.0;
                                    break;
                                case 2: n_[d] = -1.0;
                                    break;
                                default: assert ( false );
                                    break;
                            }
                        }
                        break;

                    case 5: // top
                        for ( size_t d = 0; d < DIM; ++d )
                        {
                            switch ( d )
                            {
                                case 0: n_[d] = 0.0;
                                    break;
                                case 1: n_[d] = 0.0;
                                    break;
                                case 2: n_[d] = 1.0;
                                    break;
                                default: assert ( false );
                                    break;
                            }
                        }
                        break;

                    case 1: // front
                        for ( size_t d = 0; d != DIM; ++d )
                        {
                            switch ( d )
                            {
                                case 0: n_[d] = 0.0;
                                    break;
                                case 1: n_[d] = -1.0;
                                    break;
                                case 2: n_[d] = 0.0;
                                    break;
                                default: assert ( false );
                                    break;
                            }
                        }
                        break;

                    case 4: // back
                        for ( size_t d = 0; d != DIM; ++d )
                        {
                            switch ( d )
                            {
                                case 0: n_[d] = 0.0;
                                    break;
                                case 1: n_[d] = 1.0;
                                    break;
                                case 2: n_[d] = 0.0;
                                    break;
                                default: assert ( false );
                                    break;
                            }
                        }
                        break;

                    case 2: // left
                        for ( size_t d = 0; d != DIM; ++d )
                        {
                            switch ( d )
                            {
                                case 0: n_[d] = -1.0;
                                    break;
                                case 1: n_[d] = 0.0;
                                    break;
                                case 2: n_[d] = 0.0;
                                    break;
                                default: assert ( false );
                                    break;
                            }
                        }
                        break;

                    case 3: // right
                        for ( size_t d = 0; d != DIM; ++d )
                        {
                            switch ( d )
                            {
                                case 0: n_[d] = 1.0;
                                    break;
                                case 1: n_[d] = 0.0;
                                    break;
                                case 2: n_[d] = 0.0;
                                    break;
                                default: assert ( false );
                                    break;
                            }
                        }
                        break;

                    default:
                        assert ( 0 );
                }
                break;

            case CellType::TETRAHEDRON:
                switch ( facet_number_ )
                {
                    case 0: // bottom
                        for ( size_t d = 0; d != DIM; ++d )
                        {
                            switch ( d )
                            {
                                case 0: n_[d] = 0.0;
                                    break;
                                case 1: n_[d] = 0.0;
                                    break;
                                case 2: n_[d] = -1.0;
                                    break;
                                default: assert ( false );
                                    break;
                            }
                        }
                        break;

                    case 1: // front
                        for ( size_t d = 0; d != DIM; ++d )
                        {
                            switch ( d )
                            {
                                case 0: n_[d] = 0.0;
                                    break;
                                case 1: n_[d] = -1.0;
                                    break;
                                case 2: n_[d] = 0.0;
                                    break;
                                default: assert ( false );
                                    break;
                            }
                        }
                        break;

                    case 2: // left
                        for ( size_t d = 0; d != DIM; ++d )
                        {
                            switch ( d )
                            {
                                case 0: n_[d] = -1.0;
                                    break;
                                case 1: n_[d] = 0.0;
                                    break;
                                case 2: n_[d] = 0.0;
                                    break;
                                default: assert ( false );
                                    break;
                            }
                        }
                        break;

                    case 3: // back
                        for ( size_t d = 0; d != DIM; ++d )
                        {
                            switch ( d )
                            {
                                case 0: n_[d] = 1. / std::sqrt ( 3. );
                                    break;
                                case 1: n_[d] = 1. / std::sqrt ( 3. );
                                    break;
                                case 2: n_[d] = 1. / std::sqrt ( 3. );
                                    break;
                                default: assert ( false );
                                    break;
                            }
                        }
                        break;

                    default:
                        std::cerr << "Failed: Facet number " << facet_number_ << " is not valid.\n";
                        assert ( 0 );
                }
                break;

            case CellType::PYRAMID:
                switch ( facet_number_ )
                {
                    case 0: // bottom
                        for ( size_t d = 0; d != DIM; ++d )
                        {
                            switch ( d )
                            {
                                case 0: n_[d] = 0.0;
                                    break;
                                case 1: n_[d] = 0.0;
                                    break;
                                case 2: n_[d] = -1.0;
                                    break;
                                default: assert ( false );
                                    break;
                            }
                        }
                        break;

                    case 1: // front
                        for ( size_t d = 0; d != DIM; ++d )
                        {
                            switch ( d )
                            {
                                case 0: n_[d] = 0.0;
                                    break;
                                case 1: n_[d] = -2. / std::sqrt ( 3. );
                                    break;
                                case 2: n_[d] = 1. / std::sqrt ( 3. );
                                    break;
                                default: assert ( false );
                                    break;
                            }
                        }
                        break;

                    case 2: // right
                        for ( size_t d = 0; d != DIM; ++d )
                        {
                            switch ( d )
                            {
                                case 0: n_[d] = 2. / std::sqrt ( 3. );
                                    break;
                                case 1: n_[d] = 0.0;
                                    break;
                                case 2: n_[d] = 1. / std::sqrt ( 3. );
                                    break;
                                default: assert ( false );
                                    break;
                            }
                        }
                        break;

                    case 3: // back
                        for ( size_t d = 0; d != DIM; ++d )
                        {
                            switch ( d )
                            {
                                case 0: n_[d] = 0.0;
                                    break;
                                case 1: n_[d] = 2. / std::sqrt ( 3. );
                                    break;
                                case 2: n_[d] = 1. / std::sqrt ( 3. );
                                    break;
                                default: assert ( false );
                                    break;
                            }
                        }
                        break;

                    case 4: // left
                        for ( size_t d = 0; d != DIM; ++d )
                        {
                            switch ( d )
                            {
                                case 0: n_[d] = -2. / std::sqrt ( 3. );
                                    break;
                                case 1: n_[d] = 0.0;
                                    break;
                                case 2: n_[d] = 1. / std::sqrt ( 3. );
                                    break;
                                default: assert ( false );
                                    break;
                            }
                        }
                        break;

                    default:
                        std::cerr << "Failed: Facet number " << facet_number_ << " is not valid.\n";
                        assert ( 0 );
                }
                break;

            default:
                // not yet supported
                std::cerr << "Failed: Cell type " << cell_type_ << " not yet implemented!\n";
                assert ( false );
        }

        // precompute the normals on the quadrature points.
        mapped_n_.compute ( JinvT_, EvalMappedNormal<DIM, DataType>( n_ ) );
    }

    template<int DIM, class DataType>
    template<class VectorType>
    std::vector<DataType> AssemblyAssistant<DIM, DataType>::extract_dof_values ( const VectorType& coefficients,
                                                                                 int var ) const
    {
        // extract dof-values corresponding to local variable

        // find global dof indices corresponding to this variable
        const size_t num_var_dofs = num_dofs ( var );
        std::vector<DataType> local_coefficients ( num_var_dofs, 0. );
        std::vector<int> var_dof_indices ( num_var_dofs, -1 );
        for ( size_t i = 0; i != num_var_dofs; ++i )
        {
            var_dof_indices[i] = this->global_dof_indices_[dof_index ( i, var )];
        }
        coefficients.GetValues ( vec2ptr ( var_dof_indices ), num_var_dofs, vec2ptr ( local_coefficients ) );
        return local_coefficients;
    }

}

#endif /* _ASSEMBLY_ASSISTANT_H_ */
