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

/// \author Staffan Ronnas, Simon Gawlok

#ifndef _ASSEMBLY_ASSISTANT_VALUES_H_
#    define _ASSEMBLY_ASSISTANT_VALUES_H_

#    include "assembly/function_values.h"
#    include "common/pointers.h"
#    include "common/log.h"
#    include "common/vector_algebra.h"
#    include "fem/cell_transformation.h"
#    include "fem/fetype.h"
#    include "mesh/types.h"
#    include "quadrature/quadrature.h"
#    include "space/element.h"

#    include <algorithm>
#    include <cmath>
#    include <numeric>

namespace hiflow
{
    using doffem::FEType;
    using doffem::CellTransformation;

    ////////////////////////////////////////////////////////////////////////
    //////////////// Helper functions for AssemblyAssistant ////////////////
    ////////////////////////////////////////////////////////////////////////

    // A large part of the functionality of the AssemblyAssistant is
    // implemented using the FunctionValues class evaluated with
    // different functions. The definition of these functions follows below.

    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////

    /// \brief Function that evaluates shape functions
    ///
    /// \details Evaluates all shape functions for a provided set of
    /// FEType objects. Takes as input a set of points and returns for
    /// each point, a std::vector<double> with the values of all the
    /// shape functions at that point. This is used to compute the
    /// shape function values on the reference cell.

    template<int DIM, class DataType>
    class EvalShapeFunctions
    {
      public:

        EvalShapeFunctions ( const std::vector<const FEType<DataType>* >& fe_types )
        : fe_types_ ( fe_types )
        {
        }

        inline void operator() ( int i, Vec<DIM, DataType> pt,
                std::vector<DataType>& shape_function_values ) const
        {
            shape_function_values.clear ( );
            std::vector<DataType> point ( DIM, 0.0 );
#    pragma clang loop vectorize(enable)
            for ( size_t i = 0; i < DIM; ++i )
            {
                point[i] = pt[i];
            }
            std::vector<DataType> values_for_var;

            typedef typename std::vector<const FEType<DataType>* >::const_iterator Iterator;
            for ( Iterator it = fe_types_.begin ( ), e_it = fe_types_.end ( ); it != e_it; ++it )
            {
                // no need to clear values_for_var, since it will be overwritten in N()
                values_for_var.resize ( ( *it )->get_nb_dof_on_cell ( ) );
                ( *it )->N ( point, values_for_var );
                shape_function_values.insert ( shape_function_values.end ( ),
                                               values_for_var.begin ( ),
                                               values_for_var.end ( ) );
            }
        }
      private:
        const std::vector<const FEType<DataType>* >& fe_types_;
    };

    /// \brief Function that evaluates shape function gradients
    ///
    /// \details Evaluates the gradients of all shape functions for a
    /// provided set of FEType objects. Takes as input a set of points
    /// and returns for each point, a std::vector< Vec<DIM> > with the
    /// values of all the shape function gradients at that point. This
    /// is used to compute the shape function gradients on the
    /// reference cell.

    template<int DIM, class DataType>
    class EvalShapeFunctionGradients
    {
      public:

        EvalShapeFunctionGradients ( const std::vector<const FEType<DataType>* >& fe_types )
        : fe_types_ ( fe_types )
        {
        }

        inline void operator() ( int i, Vec<DIM, DataType> pt,
                std::vector< Vec<DIM, DataType> >& shape_function_gradients ) const
        {
            shape_function_gradients.clear ( );
            std::vector<DataType> point ( DIM, 0.0 );
#    pragma clang loop vectorize(enable)
            for ( size_t i = 0; i < DIM; ++i )
            {
                point[i] = pt[i];
            }

            std::vector< Vec<DIM, DataType> > gradients_for_var;
            std::vector<DataType> gradient_component;
            typedef typename std::vector<const FEType<DataType>* >::const_iterator Iterator;
            for ( Iterator it = fe_types_.begin ( ), e_it = fe_types_.end ( ); it != e_it; ++it )
            {
                const int num_dofs_for_var = ( *it )->get_nb_dof_on_cell ( );
                gradients_for_var.clear ( );
                gradients_for_var.resize ( num_dofs_for_var );
                gradient_component.clear ( );
                gradient_component.resize ( num_dofs_for_var );

                for ( size_t c = 0; c < DIM; ++c )
                {
                    // work-around use of different functions for different gradient components
                    switch ( c )
                    {
                        case 0: ( *it )->N_x ( point, gradient_component );
                            break;
                        case 1: ( *it )->N_y ( point, gradient_component );
                            break;
                        case 2: ( *it )->N_z ( point, gradient_component );
                            break;
                        default: assert ( false );
                            break;
                    };
#    pragma clang loop vectorize(enable)
                    for ( size_t i = 0; i < num_dofs_for_var; ++i )
                    {
                        gradients_for_var[i][c] = gradient_component[i];
                    }
                }
                shape_function_gradients.insert ( shape_function_gradients.end ( ),
                                                  gradients_for_var.begin ( ),
                                                  gradients_for_var.end ( ) );
            }
        }
      private:
        const std::vector<const FEType<DataType>* >& fe_types_;
    };

    /// \brief Function that evaluates shape function hessians.
    ///
    /// \details Evaluates the hessians of all shape functions for a
    /// provided set of FEType objects. Takes as input a set of points
    /// and returns for each point, a std::vector< Mat<DIM, DIM> > with the
    /// values of all the shape function hessians at that point. This
    /// is used to compute the shape function hessians on the
    /// reference cell.

    template<int DIM, class DataType>
    class EvalShapeFunctionHessians
    {
      public:

        EvalShapeFunctionHessians ( const std::vector<const FEType<DataType>* >& fe_types )
        : fe_types_ ( fe_types )
        {
        }

        inline void operator() ( int i, Vec<DIM, DataType> pt,
                std::vector< Mat<DIM, DIM, DataType> >& shape_function_hessians ) const
        {

            shape_function_hessians.clear ( );
            std::vector<DataType> point ( DIM, 0.0 );
#    pragma clang loop vectorize(enable)
            for ( int i = 0; i < DIM; ++i )
            {
                point[i] = pt[i];
            }

            std::vector< Mat<DIM, DIM, DataType> > hessians_for_var;
            std::vector<DataType> hessian_component;
            typedef typename std::vector<const FEType<DataType>* >::const_iterator Iterator;
            for ( Iterator it = fe_types_.begin ( ), e_it = fe_types_.end ( ); it != e_it; ++it )
            {
                const int num_dofs_for_var = ( *it )->get_nb_dof_on_cell ( );
                hessians_for_var.clear ( );
                hessians_for_var.resize ( num_dofs_for_var );
                hessian_component.clear ( );
                hessian_component.resize ( num_dofs_for_var );

                for ( size_t c1 = 0; c1 < DIM; ++c1 )
                {
                    for ( size_t c2 = c1; c2 < DIM; ++c2 )
                    {
                        // work-around use of different functions for different hessian components
                        switch ( c1 )
                        {
                            case 0:
                                switch ( c2 )
                                {
                                    case 0: ( *it )->N_xx ( point, hessian_component );
                                        break;
                                    case 1: ( *it )->N_xy ( point, hessian_component );
                                        break;
                                    case 2: ( *it )->N_xz ( point, hessian_component );
                                        break;
                                    default: assert ( false );
                                        break;
                                }
                                break;
                            case 1:
                                switch ( c2 )
                                {
                                    case 1: ( *it )->N_yy ( point, hessian_component );
                                        break;
                                    case 2: ( *it )->N_yz ( point, hessian_component );
                                        break;
                                    default: assert ( false );
                                        break;
                                }
                                break;

                            case 2:
                                if ( c2 == 2 )
                                {
                                    ( *it )->N_zz ( point, hessian_component );
                                }
                                else
                                {
                                    assert ( false );
                                    break;
                                }
                                break;
                            default: assert ( false );
                                break;
                        };

                        for ( size_t i = 0; i < num_dofs_for_var; ++i )
                        {
                            hessians_for_var[i]( c1, c2 ) = hessian_component[i];
                            hessians_for_var[i]( c2, c1 ) = hessian_component[i];
                        }
                    }
                }

                shape_function_hessians.insert ( shape_function_hessians.end ( ),
                                                 hessians_for_var.begin ( ),
                                                 hessians_for_var.end ( ) );
            }
        }
      private:
        const std::vector<const FEType<DataType>* >& fe_types_;
    };

    /// \brief Function that evaluates the hessians of the cell
    /// transformation.
    ///
    /// \details Evaluates the hessians of the cell transformation at
    /// a set of points. Returns an array Mat<DIM, DIM>[DIM] with the
    /// values \f$\partial_j\partial_k\f$ of the hessian of each
    /// component \f$F_i\f$ of the cell transformation at that point.

    template<int DIM, class DataType>
    class EvalCellTransformationHessian
    {
      public:

        EvalCellTransformationHessian ( const CellTransformation<DataType>& transform )
        : transform_ ( transform )
        {
        }

        inline void operator() ( int i, Vec<DIM, DataType> pt,
                std::vector< Mat<DIM, DIM, DataType> >& hessian ) const
        {

            std::vector<DataType> point ( DIM, 0.0 );
#    pragma clang loop vectorize(enable)
            for ( size_t i = 0; i < DIM; ++i )
            {
                point[i] = pt[i];
            }
            hessian.resize ( DIM );

            hessian[0]( 0, 0 ) = transform_.x_xx ( point );

            if ( DIM == 1 ) return;

            hessian[0]( 0, 1 ) = transform_.x_xy ( point );
            hessian[0]( 1, 0 ) = hessian[0]( 0, 1 );
            hessian[0]( 1, 1 ) = transform_.x_yy ( point );

            hessian[1]( 0, 0 ) = transform_.y_xx ( point );
            hessian[1]( 0, 1 ) = transform_.y_xy ( point );
            hessian[1]( 1, 0 ) = hessian[1]( 0, 1 );
            hessian[1]( 1, 1 ) = transform_.y_yy ( point );

            if ( DIM == 2 ) return;

            hessian[0]( 0, 2 ) = transform_.x_xz ( point );
            hessian[0]( 2, 0 ) = hessian[0]( 0, 2 );
            hessian[0]( 1, 2 ) = transform_.x_yz ( point );
            hessian[0]( 2, 1 ) = hessian[0]( 1, 2 );
            hessian[0]( 2, 2 ) = transform_.x_zz ( point );

            hessian[1]( 0, 2 ) = transform_.y_xz ( point );
            hessian[1]( 2, 0 ) = hessian[1]( 0, 2 );
            hessian[1]( 1, 2 ) = transform_.y_yz ( point );
            hessian[1]( 2, 1 ) = hessian[1]( 1, 2 );
            hessian[1]( 2, 2 ) = transform_.y_zz ( point );

            hessian[2]( 0, 0 ) = transform_.z_xx ( point );
            hessian[2]( 0, 1 ) = transform_.z_xy ( point );
            hessian[2]( 1, 0 ) = hessian[2]( 0, 1 );
            hessian[2]( 1, 1 ) = transform_.z_yy ( point );

            hessian[2]( 0, 2 ) = transform_.z_xz ( point );
            hessian[2]( 2, 0 ) = hessian[2]( 0, 2 );
            hessian[2]( 1, 2 ) = transform_.z_yz ( point );
            hessian[2]( 2, 1 ) = hessian[2]( 1, 2 );
            hessian[2]( 2, 2 ) = transform_.z_zz ( point );
        }
      private:
        const CellTransformation<DataType>& transform_;
    };

    /// \brief Function that evaluates H-mapped gradients.
    ///
    /// \details Evaluates the mapping \f$H_F \nabla{}\phi\f$, where \f$H_F\f$ is
    /// the hessian of the cell transformation, and \f$\nabla{}\phi\f$ is
    /// the shape function gradient on the physical cell. Takes as
    /// input a set of points and returns for each point, a
    /// std::vector< Mat<DIM, DIM> > with the values of all the
    /// H-mapped gradients at those points. These objects are defined on the
    /// physical cell.
    ///
    /// \f$H_F\f$ is the hessian of a vector-valued function, and
    /// hence a third order tensor with components \f$H_{ijk} =
    /// \partial_j\partial_k{}H_i\f$. This function computes the
    /// reduction \f$\sum_i{H_{ijk} \cdot \partial_i\phi}\f$.

    template<int DIM, class DataType>
    class EvalHMappedGradients
    {
      public:

        EvalHMappedGradients ( const FunctionValues< std::vector< Vec<DIM, DataType> > >& grad_phi )
        : grad_phi_ ( grad_phi )
        {
        }

        inline void operator() ( int q, std::vector< Mat<DIM, DIM, DataType> > H,
                std::vector< Mat<DIM, DIM, DataType> >& h_mapped_gradients ) const
        {

            const size_t num_gradients = grad_phi_[q].size ( );
            h_mapped_gradients.resize ( num_gradients );

#    pragma clang loop vectorize(enable)
            for ( size_t i = 0; i < num_gradients; ++i )
            {
                for ( size_t c = 0; c < DIM; ++c )
                {
                    //h_mapped_gradients[i] += H[c] * grad_phi_[q][i][c];
                    h_mapped_gradients[i].Axpy ( H[c], grad_phi_[q][i][c] );
                }
            }
        }

      private:
        const FunctionValues< std::vector< Vec<DIM, DataType> > >& grad_phi_;
    };

    /// \brief Function for computing mapped shape function hessians.
    ///
    /// \details Evaluates the hessians of the shape functions on the
    /// physical element through the relation
    /// \f$H_{\phi} = J^{-T}(H_{\hat{\phi}} - H_F\nabla{\phi})J^{-1}\f$.

    template<int DIM, class DataType>
    class EvalMappedShapeFunctionHessians
    {
      public:

        EvalMappedShapeFunctionHessians (
                                          const FunctionValues< std::vector< Mat<DIM, DIM, DataType> > >& H_phi_hat,
                                          const FunctionValues< std::vector< Mat<DIM, DIM, DataType> > >& H_mapped_grad )
        : H_phi_hat_ ( H_phi_hat ), H_mapped_grad_ ( H_mapped_grad )
        {
        }

        inline void operator() ( int q, const Mat<DIM, DIM, DataType>& JinvT,
                std::vector< Mat<DIM, DIM, DataType> >& mapped_shape_function_hessians ) const
        {
            const size_t num_hessians = H_phi_hat_[q].size ( );

            mapped_shape_function_hessians.resize ( num_hessians, Mat<DIM, DIM, DataType>( ) );

            Mat<DIM, DIM, DataType> Jinv;
            trans ( JinvT, Jinv );

#    pragma clang loop vectorize(enable)
            for ( size_t i = 0; i < num_hessians; ++i )
            {
                /*mapped_shape_function_hessians[i] =
                        JinvT * (H_phi_hat_[q][i] - H_mapped_grad_[q][i]) * Jinv;*/
                // temp represents  (H_phi_hat_[q][i] - H_mapped_grad_[q][i]) * Jinv
                Mat<DIM, DIM, DataType> temp;
                for ( size_t j = 0; j < DIM; ++j )
                {
                    for ( size_t l = 0; l < DIM; ++l )
                    {
                        const DataType H_phi_temp = H_phi_hat_[q][i]( j, l ) - H_mapped_grad_[q][i]( j, l );
#    pragma clang loop vectorize(enable)
                        for ( size_t k = 0; k < DIM; ++k )
                        {
                            temp ( j, k ) += H_phi_temp * Jinv ( l, k );
                        }
                    }
                }

#    pragma clang loop vectorize(enable)
                for ( size_t j = 0; j < DIM; ++j )
                {
                    for ( size_t l = 0; l < DIM; ++l )
                    {
#    pragma clang loop vectorize(enable)
                        for ( size_t k = 0; k < DIM; ++k )
                        {
                            mapped_shape_function_hessians[i]( j, k ) += JinvT ( j, l ) * temp ( l, k );
                        }
                    }
                }
            }
        }
      private:
        const FunctionValues< std::vector< Mat<DIM, DIM, DataType> > >& H_phi_hat_;
        const FunctionValues< std::vector< Mat<DIM, DIM, DataType> > >& H_mapped_grad_;
    };

    /// \brief Function for computing mapped shape function gradients.
    ///
    /// \details Evaluates the gradients on the physical element by
    /// applying the provided set of inverse transpose:s of the
    /// jacobians of the cell transformations to each vector in the
    /// set of shape function gradients on the reference cell
    /// (grad_phi_hat). The returned vectors for each matrix JinvT are
    /// in the same order as those in grad_phi_hat.

    template<int DIM, class DataType>
    class EvalMappedShapeFunctionGradients
    {
      public:

        EvalMappedShapeFunctionGradients (
                                           const FunctionValues< std::vector< Vec<DIM, DataType> > >& grad_phi_hat )
        : grad_phi_hat_ ( grad_phi_hat )
        {
        }

        inline void operator() ( int q, const Mat<DIM, DIM, DataType>& JinvT,
                std::vector< Vec<DIM, DataType> >& mapped_shape_function_gradients ) const
        {
            const size_t num_gradients = grad_phi_hat_[q].size ( );
            mapped_shape_function_gradients.resize ( num_gradients );
            for ( size_t i = 0; i < num_gradients; ++i )
            {
                //mapped_shape_function_gradients[i] = JinvT * grad_phi_hat_[q][i];
                JinvT.VectorMult ( grad_phi_hat_[q][i], mapped_shape_function_gradients[i] );
            }
        }
      private:
        const FunctionValues< std::vector< Vec<DIM, DataType> > >& grad_phi_hat_;
    };

    /// \brief Maps a normal of a reference element to the
    /// physical element.

    template<int DIM, class DataType>
    class EvalMappedNormal
    {
      public:

        EvalMappedNormal (
                           const Vec<DIM, DataType>& n )
        : n_ ( n )
        {
        }

        inline void operator() ( int q, const Mat<DIM, DIM, DataType>& JinvT,
                Vec<DIM, DataType>& mapped_n ) const
        {
            //mapped_n = JinvT * n_;
            JinvT.VectorMult ( n_, mapped_n );
            mapped_n /= norm ( mapped_n );
        }
      private:
        const Vec<DIM, DataType>& n_;
    };
    /// \brief Computes the coordinates of a set of points on the physical element.
    ///
    /// \details Transforms a set of points on the reference element to
    /// the physical element, using the provided CellTransform
    /// object. This is used to compute the coordinates of the
    /// quadrature points on the physical element.

    template<int DIM, class DataType>
    class EvalPhysicalPoint
    {
      public:

        EvalPhysicalPoint ( const CellTransformation<DataType>& transform )
        : transform_ ( transform )
        {
            assert ( DIM > 0 );
            assert ( DIM <= 3 );
        }

        inline void operator() ( int i, Vec<DIM, DataType> pt, Vec<DIM, DataType>& mapped_point ) const
        {
            std::vector<DataType> point ( DIM, 0.0 );
#    pragma clang loop vectorize(enable)
            for ( size_t i = 0; i < DIM; ++i )
            {
                point[i] = pt[i];
            }

            switch ( DIM )
            {
                case 0:
                {
                    assert ( false );
                    break;
                }
                case 1:
                {
                    mapped_point[0] = transform_.x ( point );
                    break;
                }
                case 2:
                {
                    mapped_point[0] = transform_.x ( point );
                    mapped_point[1] = transform_.y ( point );
                    break;
                }
                case 3:
                {
                    mapped_point[0] = transform_.x ( point );
                    mapped_point[1] = transform_.y ( point );
                    mapped_point[2] = transform_.z ( point );
                    break;
                }
                default:
                {
                    assert ( false );
                    break;
                }
            }
        }
      private:
        const CellTransformation<DataType>& transform_;
    };

    /// \brief Computes the jacobian of the element transformation at
    /// a set of points on the reference element.
    ///
    /// \details This is used to compute the jacobian of the element
    /// transformation at the quadrature points.

    template<int DIM, class DataType>
    class EvalPhysicalJacobian
    {
      public:

        EvalPhysicalJacobian ( const CellTransformation<DataType>& transform )
        : transform_ ( transform )
        {
            assert ( DIM > 0 );
            assert ( DIM <= 3 );
        }

        inline void operator() ( int i, Vec<DIM, DataType> pt, Mat<DIM, DIM, DataType>& jacobian ) const
        {
            std::vector<DataType> point ( DIM, 0.0 );
#    pragma clang loop vectorize(enable)
            for ( size_t i = 0; i < DIM; ++i )
            {
                point[i] = pt[i];
            }

            jacobian ( 0, 0 ) = transform_.x_x ( point );
            if ( DIM > 1 )
            {
                jacobian ( 0, 1 ) = transform_.x_y ( point );
                jacobian ( 1, 0 ) = transform_.y_x ( point );
                jacobian ( 1, 1 ) = transform_.y_y ( point );
            }

            if ( DIM > 2 )
            {
                jacobian ( 0, 2 ) = transform_.x_z ( point );
                jacobian ( 1, 2 ) = transform_.y_z ( point );

                jacobian ( 2, 0 ) = transform_.z_x ( point );
                jacobian ( 2, 1 ) = transform_.z_y ( point );
                jacobian ( 2, 2 ) = transform_.z_z ( point );
            }

        }
      private:
        const CellTransformation<DataType>& transform_;
    };

    /// \brief Computes the determinant of a set of matrices.
    ///
    /// \details Each evaluation computes the determinant of a
    /// matrix. This is used to compute the determinants of the
    /// jacobian matrices at the quadrature points.

    template<int DIM, class DataType>
    struct EvalDeterminant
    {

        inline void operator() ( int i, const Mat<DIM, DIM, DataType>& matrix, DataType& determinant ) const
        {
            determinant = det ( matrix );
        }
    };

    /// \brief Computes the inverse transpose of a set of matrices.
    ///
    /// \details Each evaluation computes the inverse transpose of a
    /// matrix. This is used to compute the inverse transposes of the
    /// jacobian matrices at the quadrature points.

    template<int DIM, class DataType>
    struct EvalInvTranspose
    {

        void operator() ( int i, const Mat<DIM, DIM, DataType>& matrix, Mat<DIM, DIM, DataType>& inv_T ) const
        {
            invTransp ( matrix, inv_T );
        }
    };

    /// \brief Multiplies matrices on the right with a given matrix.
    ///
    /// \details Each evaluation computes B = A*R. The dimensions of
    /// the matrices are as follows: A -> MxN, R -> NxP, B->MxP.

    template<int M, int N, int P, class DataType>
    struct EvalRightMatrixMult
    {

        EvalRightMatrixMult ( const Mat<N, P, DataType>& R )
        : R_ ( R )
        {
        }

        inline void operator() ( int i, const Mat<M, N, DataType>& A, Mat<M, P, DataType>& B ) const
        {
            //B = A * R_;
            MatrixMatrixMult ( B, A, R_ );
        }

      private:
        const Mat<N, P, DataType>& R_;
    };

    /// \brief Evaluate surface element.
    ///
    /// \brief Given Jacobian matrix Jf of mapping R^{D-1} -> R^D, computes
    /// the surface element ds = \sqrt(det(Jf^T * Jf)).

    template<int DIM, class DataType>
    struct EvalSurfaceElement
    {

        inline void operator() ( int i, const Mat<DIM, DIM - 1, DataType>& Jf, DataType& ds ) const
        {
            Mat < DIM - 1, DIM, DataType> JfT;
            trans ( Jf, JfT );

            Mat < DIM - 1, DIM - 1, DataType> JfTJf;
            MatrixMatrixMult ( JfTJf, JfT, Jf );

            const DataType detJfTJf = det ( JfTJf );
            assert ( detJfTJf > 0. );

            ds = std::sqrt ( detJfTJf );
        }
    };

    /// \brief Evaluates a finite element function defined through the
    /// values of its degrees of freedoms for different sets of shape
    /// function values, typically corresponding to different points.
    ///
    /// \details EvalFiniteElementFunction takes as input a set of
    /// shape function values {\phi_i}, and evaluates \sum{u_i *
    /// \phi_i}. The values of u_i are given in the variable
    /// local_coefficients. The index set of i is offset by the
    /// variable fe_offset, so that i goes from
    ///
    /// [fe_offset, ..., fe_offset + local_coefficients.size() [ .
    ///
    /// The offset makes it possible to compute the values only for an
    /// isolated variable. This function is used in the context of
    /// AssemblyAssistant::evaluate_fe_function().

    template<class DataType>
    class EvalFiniteElementFunction
    {
      public:

        EvalFiniteElementFunction ( const int fe_offset, const std::vector<DataType>& local_coefficients )
        : fe_offset_ ( fe_offset ), local_coefficients_ ( local_coefficients )
        {
        }

        inline void operator() ( int q, const std::vector<DataType>& phi_values, DataType& u_value ) const
        {
            const size_t num_dofs = local_coefficients_.size ( );

            u_value = 0.;
#    pragma clang loop vectorize(enable)
            for ( size_t i = 0; i < num_dofs; ++i )
            {
                u_value += local_coefficients_[i] * phi_values[fe_offset_ + i];
            }
        }
      private:
        const int fe_offset_;
        const std::vector<DataType>& local_coefficients_;
    };

    /// \brief Evaluates the gradient of a finite element function
    /// defined through the values of its degrees of freedoms for
    /// different sets of shape function values, typically
    /// corresponding to different points.
    ///
    /// \details Does the same as \see EvalFiniteElementFunction in
    /// principle, but computes \sum{u_i * \nabla{phi}_i} instead. This
    /// is used in the context of
    /// AssemblyAssistant::evaluate_fe_function_gradients().

    template<int DIM, class DataType>
    class EvalFiniteElementFunctionGradient
    {
      public:

        EvalFiniteElementFunctionGradient ( const int fe_offset,
                                            const std::vector<DataType>& local_coefficients )
        : fe_offset_ ( fe_offset ), local_coefficients_ ( local_coefficients )
        {
        }

        inline void operator() ( int q, const std::vector< Vec<DIM, DataType> >& grad_phi_values,
                Vec<DIM, DataType>& grad_u_value ) const
        {
            const size_t num_dofs = local_coefficients_.size ( );

            grad_u_value = Vec<DIM, DataType>( );
            for ( size_t i = 0; i != num_dofs; ++i )
            {
                //grad_u_value += local_coefficients_[i] * grad_phi_values[fe_offset_ + i];
                grad_u_value.Axpy ( grad_phi_values[fe_offset_ + i], local_coefficients_[i] );
            }
        }
      private:
        const int fe_offset_;
        const std::vector<DataType>& local_coefficients_;
    };

    /// \brief Evaluates the hessian of a finite element function
    /// defined through the values of its degrees of freedoms for
    /// different sets of shape function values, typically
    /// corresponding to different points.
    ///
    /// \details Does the same as \see EvalFiniteElementFunction in
    /// principle, but computes \sum{u_i * H(phi_i)} instead. This
    /// is used in the context of
    /// AssemblyAssistant::evaluate_fe_function_hessians().

    template<int DIM, class DataType>
    class EvalFiniteElementFunctionHessian
    {
      public:

        EvalFiniteElementFunctionHessian ( const int fe_offset,
                                           const std::vector<DataType>& local_coefficients )
        : fe_offset_ ( fe_offset ), local_coefficients_ ( local_coefficients )
        {
        }

        inline void operator() ( int q, const std::vector< Mat<DIM, DIM, DataType> >& H_phi_values,
                Mat<DIM, DIM, DataType>& H_u_value ) const
        {
            const size_t num_dofs = local_coefficients_.size ( );

            H_u_value = Mat<DIM, DIM, DataType>( );
            for ( size_t i = 0; i != num_dofs; ++i )
            {
                //H_u_value += local_coefficients_[i] * H_phi_values[fe_offset_ + i];
                H_u_value.Axpy ( H_phi_values[fe_offset_ + i], local_coefficients_[i] );
            }
        }
      private:
        const int fe_offset_;
        const std::vector<DataType>& local_coefficients_;
    };

}

#endif /* _ASSEMBLY_ASSISTANT_VALUES_H_ */
