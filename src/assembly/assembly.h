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

#ifndef _ASSEMBLY_H_
#    define _ASSEMBLY_H_

#    include "linear_algebra/vector.h"
#    include "linear_algebra/matrix.h"
#    include "linear_algebra/couplings.h"
#    include "linear_algebra/seq_dense_matrix.h"
#    include "quadrature/quadrature.h"
#    include "space/element.h"
#    include "space/solution.h"
#    include "space/vector_space.h"

#    include <fstream>
#    include <vector>

#    include <boost/function.hpp>

/// \file assembly.h
/// \brief Assembly functions.
///
/// \author Staffan Ronnas
///

namespace hiflow
{

    namespace la
    {
        template<class DataType>
        class Couplings;
    }

    template<class DataType>
    class Element;

    template<class DataType>
    class VectorSpace;

    struct SparsityStructure
    {
        std::vector<int> diagonal_rows;
        std::vector<int> diagonal_cols;
        std::vector<int> off_diagonal_rows;
        std::vector<int> off_diagonal_cols;
    };

    /// Abstract interface for classes implementing a Global Assembly strategy.

    template<class DataType>
    class GlobalAssembler
    {
      public:
        typedef la::Matrix<DataType> GlobalMatrix;
        typedef la::Vector<DataType> GlobalVector;
        typedef la::SeqDenseMatrix<DataType> LocalMatrix;
        typedef std::vector<DataType> LocalVector;

        /// Generalized function type used for scalar local assembly functions.
        typedef boost::function3<void,
        const Element<DataType>&,
        const Quadrature<DataType>&,
        DataType&> ScalarAssemblyFunction;

        /// Generalized function type used for multiple scalar local assembly functions.
        typedef boost::function3<void,
        const Element<DataType>&,
        const Quadrature<DataType>&,
        LocalVector&> MultipleScalarAssemblyFunction;

        /// Generalized function type used for vector local assembly functions.
        typedef boost::function3<void,
        const Element<DataType>&,
        const Quadrature<DataType>&,
        LocalVector&> VectorAssemblyFunction;

        /// Generalized function type used for matrix local assembly functions.
        typedef boost::function3<void,
        const Element<DataType>&,
        const Quadrature<DataType>&,
        LocalMatrix&> MatrixAssemblyFunction;

        /// Generalized function type used for scalar local assembly functions on the boundary.
        typedef boost::function4<void,
        const Element<DataType>&,
        int,
        const Quadrature<DataType>&,
        LocalVector&> BoundaryScalarAssemblyFunction;

        /// Generalized function type used for vector local assembly functions on the boundary.
        typedef boost::function4<void,
        const Element<DataType>&,
        int,
        const Quadrature<DataType>&,
        LocalVector&> BoundaryVectorAssemblyFunction;

        /// Generalized function type used for matrix local assembly functions on the boundary.
        typedef boost::function4<void,
        const Element<DataType>&,
        int,
        const Quadrature<DataType>&,
        LocalMatrix&> BoundaryMatrixAssemblyFunction;

        /// Generalized function type used for quadrature selection.
        typedef boost::function2<void,
        const Element<DataType>&,
        Quadrature<DataType>&> QuadratureSelectionFunction;

        typedef boost::function3<void,
        const Element<DataType>&,
        Quadrature<DataType>&, const int> FacetQuadratureSelectionFunction;

        GlobalAssembler ( );
        virtual ~GlobalAssembler ( );

        ///
        /// \brief Compute the matrix graph for the assembly.
        ///
        /// \details The matrix graph is a set of pairs (i,j) of
        /// global indices which are cannot be guaranteed to give a
        /// zero in the matrix assembly. Typically, these will
        /// correspond to all global basis functions with overlapping
        /// support.
        ///
        /// \param[in]  space        the VectorSpace to assemble over
        /// \param[out] sparsity     the sparsity object containing
        ///                          arrays that describe the matrix graph
        /// \param[in] coupling_vars 2D array indicating the coupling of vars.
        ///                          Rows (first index) belong to test variables,
        ///                          columns (second index) belong to trial variables.
        ///                          If entry (i, j) is set to true, test variable i
        ///                          couples with trial variable j, and the corresponding
        ///                          block is contained in the sparsity structure. Otherwise,
        ///                          the block is skipped, resulting in a sparser structure.
        ///                          If this argument is not passed or is empty, full coupling
        ///                          of all variables is assumed. All rows and columns need
        ///                          to have the size of space.get_nb_var().
        void compute_sparsity_structure ( const VectorSpace<DataType>& space,
                                          SparsityStructure& sparsity,
                                          std::vector<std::vector<bool> > *coupling_vars = new std::vector<std::vector<bool> >( 0 ) ) const;

        ///
        /// \brief Compute integral for a VectorSpace.
        ///
        /// \details Computes an integral \f$\int_{\Omega}{f(x)dx}\f$
        /// over the domain defined by the mesh associated to a
        /// VectorSpace. The integrand is defined through the
        /// ScalarAssemblyFunction, which should return the value of
        /// the integral for each element.
        ///
        /// \param[in] space      the VectorSpace to integrate over
        /// \param[in] local_asm  function or functor that performs local integration
        /// \param[out] value     the value of the integral
        /// \see concept_assembly
        ///
        void integrate_scalar ( const VectorSpace<DataType>& space,
                                ScalarAssemblyFunction local_asm,
                                DataType& integral ) const;

        ///
        /// \brief Compute integral for a VectorSpace.
        ///
        /// \details Computes an integral \f$\int_{\Omega}{f(x)dx}\f$
        /// over the domain defined by the mesh associated to a
        /// VectorSpace. The vector valued integrand is defined through the
        /// MultipleScalarAssemblyFunction, which should return the value of
        /// the integral for each element.
        ///
        /// \param[in] space      the VectorSpace to integrate over
        /// \param[in] local_asm  function or functor that performs local integration
        /// \param[in] num_scalars dimension of vector valued integrand						  
        /// \param[out] value     the value of the integral
        /// \see concept_assembly
        ///
        void integrate_multiple_scalar ( const VectorSpace<DataType>& space,
                                         MultipleScalarAssemblyFunction local_asm,
                                         const int num_scalars,
                                         std::vector<DataType>& integral ) const;

        ///
        /// \brief Compute element-wise integral for a VectorSpace.
        ///
        /// \details Computes an integral \f$\int_{K}{f(x)dx}\f$
        /// over the cells \f$K\f$ in the mesh associated to a
        /// VectorSpace. The integrand is defined through the
        /// ScalarAssemblyFunction, which should return the value of
        /// the integral for each element.
        ///
        /// \param[in] space      the VectorSpace to integrate over
        /// \param[in] local_asm  function or functor that performs local integration
        /// \param[out] vec       the value of the integral for each cell in the mesh
        /// \see concept_assembly
        ///
        void assemble_scalar ( const VectorSpace<DataType>& space,
                               ScalarAssemblyFunction local_asm,
                               std::vector<DataType>& vec ) const;

        ///
        /// \brief Compute element-wise integral for a VectorSpace.
        ///
        /// \details Computes an integral \f$\int_{K}{f(x)dx}\f$
        /// over the cells \f$K\f$ in the mesh associated to a
        /// VectorSpace. The vector-valued integrand is defined through the
        /// MultipleScalarAssemblyFunction, which should return the values of
        /// the integral for each element.
        ///
        /// \param[in] space      the VectorSpace to integrate over
        /// \param[in] local_asm  function or functor that performs local integration
        /// \param[in] num_scalars dimension of vector valued integrand					 
        /// \param[out] vec       the values of the integral for each cell in the mesh
        /// \see concept_assembly
        ///
        void assemble_multiple_scalar ( const VectorSpace<DataType>& space,
                                        MultipleScalarAssemblyFunction local_asm,
                                        const int num_scalars,
                                        std::vector< std::vector<DataType> >& vec ) const;

        ///
        /// \brief Assemble global vector for a VectorSpace.
        ///
        /// \details Assembles a vector
        /// \f$b_i = \int_{\Omega}{f(x, \varphi_i)dx}\f$
        /// over the domain defined by the mesh
        /// associated to a VectorSpace. The integrand is defined through
        /// a VectorAssemblyFunction, which should return the locally
        /// assembled vector for each element.
        ///
        /// \param[in] space      the VectorSpace for which the assembly is performed
        /// \param[in] local_asm  function or functor that performs local vector assembly
        /// \param[out] global_vector  the assembled vector \f$b_i\f$
        /// \see concept_assembly
        ///
        void assemble_vector ( const VectorSpace<DataType>& space,
                               VectorAssemblyFunction local_asm,
                               GlobalVector& vec ) const;

        ///
        /// \brief Assemble global matrix for a VectorSpace.
        ///
        /// \details Assembles a matrix
        /// \f$A_{ij} = \int_{\Omega}{f(x,\varphi_i, \varphi_j)dx}\f$
        /// over the domain defined by the mesh
        /// associated to a VectorSpace. The integrand is defined through
        /// a MatrixAssemblyFunction, which should return the locally
        /// assembled matrix for each element.
        ///
        /// \param[in] space           the VectorSpace for which the assembly is performed
        /// \param[in] local_asm       function or functor that performs local matrix assembly
        /// \param[out] global_matrix  the assembled matrix \f$A_{ij}\f$
        /// \see concept_assembly
        ///
        void assemble_matrix ( const VectorSpace<DataType>& space,
                               MatrixAssemblyFunction local_asm,
                               GlobalMatrix& matrix ) const;

        ///
        /// \brief Compute boundary integral for a VectorSpace.
        ///
        /// \details Computes a boundary integral
        /// \f$\int_{\partial\Omega}{f(x)dx}\f$ over the boundary of
        /// the domain defined by the mesh associated to a
        /// VectorSpace. The integrand is defined through the
        /// BoundaryScalarAssemblyFunction, which should return the value of
        /// the integral for each element.
        ///
        /// \param[in] space      the VectorSpace to integrate over
        /// \param[in] local_asm  function or functor that performs local integration
        /// \param[out] value     the value of the integral
        /// \see concept_assembly
        ///
        void integrate_scalar_boundary ( const VectorSpace<DataType>& space,
                                         BoundaryScalarAssemblyFunction local_asm,
                                         DataType& integral ) const;

        ///
        /// \brief Compute boundary maximum for a VectorSpace.
        ///
        /// \details Determines a maximum value
        /// \f$\max_{\partial\Omega}{f(x)}\f$ at the boundary of
        /// the domain defined by the mesh associated to a
        /// VectorSpace. The cell maxima are defined through the
        /// BoundaryScalarAssemblyFunction, which should return the maximal value
        /// for each element.
        ///
        /// \param[in] space      the VectorSpace to maximize over
        /// \param[in] local_asm  function or functor that performs local maximization
        /// \param[out] maximum   the maximal value
        /// \see concept_assembly
        ///
        void maximize_scalar_boundary ( const VectorSpace<DataType>& space,
                                        BoundaryScalarAssemblyFunction local_asm,
                                        DataType& maximum ) const;

        ///
        /// \brief Compute facet-wise boundary integral for a VectorSpace.
        ///
        /// \details Computes a boundary integral \f$\int_{\partial K}{f(x)dx}\f$
        /// over the cells in the mesh associated to a
        /// VectorSpace. The integrand is defined through the
        /// BoundaryScalarAssemblyFunction, which should return the value of
        /// the integral for each boundary facet.
        ///
        /// \param[in] space      the VectorSpace to integrate over
        /// \param[in] local_asm  function or functor that performs local integration
        /// \param[out] vec       the value of the integral for each facet on the boundary
        /// \see concept_assembly
        ///
        void assemble_scalar_boundary ( const VectorSpace<DataType>& space,
                                        BoundaryScalarAssemblyFunction local_asm,
                                        std::vector<DataType>& vec ) const;

        ///
        /// \brief Assemble global vector for a VectorSpace, defined
        /// by a boundary integral.
        ///
        /// \details Assembles a vector
        /// \f$b_i = \int_{\partial\Omega}{f(x, \varphi_i)ds}\f$
        /// over the boundary of the domain defined by the mesh
        /// associated to a VectorSpace. The integrand is defined through
        /// a BoundaryVectorAssemblyFunction, which should return the locally
        /// assembled vector for each boundary facet.
        ///
        /// \param[in] space           the VectorSpace for which the assembly is performed
        /// \param[in] local_asm       function or functor that performs local vector assembly
        /// \param[out] global_vector  the assembled vector \f$b_i\f$
        /// \see concept_assembly
        ///
        void assemble_vector_boundary ( const VectorSpace<DataType>& space,
                                        BoundaryVectorAssemblyFunction local_asm,
                                        GlobalVector& vec ) const;

        ///
        /// \brief Assemble global matrix for a VectorSpace, defined
        /// by a boundary integral.
        ///
        /// \details Assembles a matrix
        /// \f$A_{ij} = \int_{\partial\Omega}{f(x,\varphi_i, \varphi_j)dx}\f$
        /// over the boundary of the domain defined by the mesh
        /// associated to a VectorSpace. The integrand is defined through
        /// a BoundaryMatrixAssemblyFunction, which should return the locally
        /// assembled matrix for each boundary facet.
        ///
        /// \param[in] space           the VectorSpace for which the assembly is performed
        /// \param[in] local_asm       function or functor that performs local matrix assembly
        /// \param[out] global_matrix  the assembled matrix \f$A_{ij}\f$
        /// \see concept_assembly
        ///
        void assemble_matrix_boundary ( const VectorSpace<DataType>& space,
                                        BoundaryMatrixAssemblyFunction local_asm,
                                        GlobalMatrix& matrix ) const;

        ///
        /// \brief Set the function used to determine which quadrature rule should be used.
        ///
        /// \details The choice of quadrature rules can be controlled
        /// by providing a QuadratureSelectionFunction. This function
        /// or functor will be called before the local assembly is
        /// performed on each element. By default, the
        /// DefaultQuadratureSelection function is used.
        ///
        /// \param[in] q_select   new quadrature selection function
        void set_quadrature_selection_function ( QuadratureSelectionFunction q_select );

        ///
        /// \brief Change whether the assembled object should be
        /// reset prior to assembly.
        ///
        /// \details This function can be used to decide whether the
        /// target object of the assembly (scalar for scalar assembly,
        /// vector for vector assembly, etc) in the various assembly
        /// functions should be reset to zero before the
        /// assembly. Setting this value to false can be used to
        /// perform several assembly steps with the same object,
        /// e.g. for adding a boundary term to an assembly. By
        /// default, this is set to true.
        ///
        /// \param[in] should_reset   whether or not assembly object should be reset to zero.
        void should_reset_assembly_target ( bool should_reset );

      private:

        /// Interface to be implemented by concrete assembly
        /// classes. These are called from the corresponding functions
        /// in the public interface. The integrate_scalar_* functions
        /// call assemble_scalar_impl and sum the values in the
        /// vector. Note that the implementations should not reset the
        /// global object themselves, as this is handled by the public
        /// functions.
        virtual void compute_sparsity_structure_impl ( const VectorSpace<DataType>& space,
                                                       SparsityStructure& sparsity,
                                                       std::vector<std::vector<bool> > *coupling_vars ) const = 0;

        virtual void assemble_scalar_impl ( const VectorSpace<DataType>& space,
                                            ScalarAssemblyFunction local_asm,
                                            std::vector<DataType>& vec,
                                            QuadratureSelectionFunction q_select ) const = 0;

        virtual void assemble_multiple_scalar_impl ( const VectorSpace<DataType>& space,
                                                     MultipleScalarAssemblyFunction local_asm,
                                                     const int num_scalars,
                                                     std::vector< std::vector<DataType> >& vec,
                                                     QuadratureSelectionFunction q_select ) const = 0;

        virtual void assemble_vector_impl ( const VectorSpace<DataType>& space,
                                            VectorAssemblyFunction local_asm,
                                            GlobalVector& vec,
                                            QuadratureSelectionFunction q_select ) const = 0;

        virtual void assemble_matrix_impl ( const VectorSpace<DataType>& space,
                                            MatrixAssemblyFunction local_asm,
                                            GlobalMatrix& mat,
                                            QuadratureSelectionFunction q_select ) const = 0;

        virtual void assemble_scalar_boundary_impl ( const VectorSpace<DataType>& space,
                                                     BoundaryScalarAssemblyFunction local_asm,
                                                     std::vector<DataType>& vec,
                                                     FacetQuadratureSelectionFunction fq_select ) const = 0;

        virtual void assemble_vector_boundary_impl ( const VectorSpace<DataType>& space,
                                                     BoundaryVectorAssemblyFunction local_asm,
                                                     GlobalVector& vec,
                                                     FacetQuadratureSelectionFunction fq_select ) const = 0;

        virtual void assemble_matrix_boundary_impl ( const VectorSpace<DataType>& space,
                                                     BoundaryMatrixAssemblyFunction local_asm,
                                                     GlobalMatrix& mat,
                                                     FacetQuadratureSelectionFunction fq_select ) const = 0;

        static const QuadratureSelectionFunction default_select_;
        static const FacetQuadratureSelectionFunction default_facet_select_;
        QuadratureSelectionFunction q_select_;
        FacetQuadratureSelectionFunction fq_select_;

      protected:
        bool should_reset_assembly_target_;
    };

    /// Standard Global Assembly strategy

    template<class DataType>
    class StandardGlobalAssembler : public GlobalAssembler<DataType>
    {
        virtual void compute_sparsity_structure_impl ( const VectorSpace<DataType>& space,
                                                       SparsityStructure& sparsity,
                                                       std::vector<std::vector<bool> > *coupling_vars ) const;

        virtual void assemble_scalar_impl ( const VectorSpace<DataType>& space,
                                            typename GlobalAssembler<DataType>::ScalarAssemblyFunction local_asm,
                                            std::vector<DataType>& vec,
                                            typename GlobalAssembler<DataType>::QuadratureSelectionFunction q_select ) const;

        virtual void assemble_multiple_scalar_impl ( const VectorSpace<DataType>& space,
                                                     typename GlobalAssembler<DataType>::MultipleScalarAssemblyFunction local_asm,
                                                     const int num_scalars,
                                                     std::vector< std::vector<DataType> >& vec,
                                                     typename GlobalAssembler<DataType>::QuadratureSelectionFunction q_select ) const;

        virtual void assemble_vector_impl ( const VectorSpace<DataType>& space,
                                            typename GlobalAssembler<DataType>::VectorAssemblyFunction local_asm,
                                            typename GlobalAssembler<DataType>::GlobalVector& vec,
                                            typename GlobalAssembler<DataType>::QuadratureSelectionFunction q_select ) const;

        virtual void assemble_matrix_impl ( const VectorSpace<DataType>& space,
                                            typename GlobalAssembler<DataType>::MatrixAssemblyFunction local_asm,
                                            typename GlobalAssembler<DataType>::GlobalMatrix& mat,
                                            typename GlobalAssembler<DataType>::QuadratureSelectionFunction q_select ) const;

        virtual void assemble_scalar_boundary_impl ( const VectorSpace<DataType>& space,
                                                     typename GlobalAssembler<DataType>::BoundaryScalarAssemblyFunction local_asm,
                                                     std::vector<DataType>& vec,
                                                     typename GlobalAssembler<DataType>::FacetQuadratureSelectionFunction fq_select ) const;

        virtual void assemble_vector_boundary_impl ( const VectorSpace<DataType>& space,
                                                     typename GlobalAssembler<DataType>::BoundaryVectorAssemblyFunction local_asm,
                                                     typename GlobalAssembler<DataType>::GlobalVector& vec,
                                                     typename GlobalAssembler<DataType>::FacetQuadratureSelectionFunction fq_select ) const;

        virtual void assemble_matrix_boundary_impl ( const VectorSpace<DataType>& space,
                                                     typename GlobalAssembler<DataType>::BoundaryMatrixAssemblyFunction local_asm,
                                                     typename GlobalAssembler<DataType>::GlobalMatrix& mat,
                                                     typename GlobalAssembler<DataType>::FacetQuadratureSelectionFunction fq_select ) const;
    };

    /// hp-FEM Global Assembly strategy

    template<class DataType>
    class HpFemAssembler : public GlobalAssembler<DataType>
    {
        virtual void compute_sparsity_structure_impl ( const VectorSpace<DataType>& space,
                                                       SparsityStructure& sparsity,
                                                       std::vector<std::vector<bool> > *coupling_vars ) const;

        virtual void assemble_scalar_impl ( const VectorSpace<DataType>& space,
                                            typename GlobalAssembler<DataType>::ScalarAssemblyFunction local_asm,
                                            std::vector<DataType>& vec,
                                            typename GlobalAssembler<DataType>::QuadratureSelectionFunction q_select ) const;

        virtual void assemble_multiple_scalar_impl ( const VectorSpace<DataType>& space,
                                                     typename GlobalAssembler<DataType>::MultipleScalarAssemblyFunction local_asm,
                                                     const int num_scalars,
                                                     std::vector< std::vector<DataType> >& vec,
                                                     typename GlobalAssembler<DataType>::QuadratureSelectionFunction q_select ) const;

        virtual void assemble_vector_impl ( const VectorSpace<DataType>& space,
                                            typename GlobalAssembler<DataType>::VectorAssemblyFunction local_asm,
                                            typename GlobalAssembler<DataType>::GlobalVector& vec,
                                            typename GlobalAssembler<DataType>::QuadratureSelectionFunction q_select ) const;

        virtual void assemble_matrix_impl ( const VectorSpace<DataType>& space,
                                            typename GlobalAssembler<DataType>::MatrixAssemblyFunction local_asm,
                                            typename GlobalAssembler<DataType>::GlobalMatrix& mat,
                                            typename GlobalAssembler<DataType>::QuadratureSelectionFunction q_select ) const;

        virtual void assemble_scalar_boundary_impl ( const VectorSpace<DataType>& space,
                                                     typename GlobalAssembler<DataType>::BoundaryScalarAssemblyFunction local_asm,
                                                     std::vector<DataType>& vec,
                                                     typename GlobalAssembler<DataType>::FacetQuadratureSelectionFunction fq_select ) const;

        virtual void assemble_vector_boundary_impl ( const VectorSpace<DataType>& space,
                                                     typename GlobalAssembler<DataType>::BoundaryVectorAssemblyFunction local_asm,
                                                     typename GlobalAssembler<DataType>::GlobalVector& vec,
                                                     typename GlobalAssembler<DataType>::FacetQuadratureSelectionFunction fq_select ) const;

        virtual void assemble_matrix_boundary_impl ( const VectorSpace<DataType>& space,
                                                     typename GlobalAssembler<DataType>::BoundaryMatrixAssemblyFunction local_asm,
                                                     typename GlobalAssembler<DataType>::GlobalMatrix& mat,
                                                     typename GlobalAssembler<DataType>::FacetQuadratureSelectionFunction fq_select ) const;
    };

    //////////////////////////////////////////////////////////////////////////////////
    //////////////// Various functions ///////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    /// \brief Initialize linear algebra objects for global assembly.
    template<class DataType>
    void initialize_linear_algebra_objects ( const VectorSpace<DataType>& space,
                                             const GlobalAssembler<DataType>& global_asm,
                                             la::Couplings<DataType>& couplings,
                                             typename GlobalAssembler<DataType>::GlobalMatrix& matrix );

    template<class DataType>
    void interpolate_constrained_vector ( const VectorSpace<DataType>& space, typename GlobalAssembler<DataType>::GlobalVector& vector );
}

#endif /* _ASSEMBLY_H_ */
