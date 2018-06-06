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

#ifndef HIFLOW_ASSEMBLY_DG_ASSEMBLY_H
#    define HIFLOW_ASSEMBLY_DG_ASSEMBLY_H

#    include "assembly/assembly.h"
#    include "linear_algebra/coupled_matrix.h"
#    include "mesh/interface.h"
#    include "quadrature/quadrature.h"

#    include <boost/function.hpp>

/// \author Staffan Ronnas, Jonathan Schwegler, Simon Gawlok

namespace hiflow
{

    /// \author Staffan Ronnas, Jonathan Schwegler, Simon Gawlok

    /// \brief Class for performing global assembly over interfaces for e.g.
    /// Discontinuous Galerkin-Methods.

    template <class DataType>
    class DGGlobalAssembler : public StandardGlobalAssembler<DataType>
    {
      public:

        enum InterfaceSide
        {
            INTERFACE_MASTER = 0, INTERFACE_SLAVE = 1, INTERFACE_BOUNDARY = 2
        };
        typedef boost::function9<void,
        const Element<DataType>&, const Element<DataType>&,
        const Quadrature<DataType>&,
        const Quadrature<DataType>&,
        int, int,
        InterfaceSide, InterfaceSide,
        typename GlobalAssembler<DataType>::LocalMatrix&
        > InterfaceMatrixAssemblyFun;

        typedef boost::function9<void,
        const Element<DataType>&, const Element<DataType>&,
        const Quadrature<DataType>&,
        const Quadrature<DataType>&,
        int, int,
        InterfaceSide, InterfaceSide,
        typename GlobalAssembler<DataType>::LocalVector&
        > InterfaceVectorAssemblyFun;

        typedef boost::function9<void,
        const Element<DataType>&, const Element<DataType>&,
        const Quadrature<DataType>&,
        const Quadrature<DataType>&,
        int, int,
        InterfaceSide, InterfaceSide,
        DataType&
        > InterfaceScalarAssemblyFun;

        typedef boost::function9<void,
        const Element<DataType>&, const Element<DataType>&,
        const Quadrature<DataType>&,
        const Quadrature<DataType>&,
        int, int,
        InterfaceSide, InterfaceSide,
        typename GlobalAssembler<DataType>::LocalVector&
        > InterfaceMultipleScalarAssemblyFun;

        // TODO: this could be optimized by splitting into separate cell/facet
        // selection functions.
        typedef boost::function6<void,
        const Element<DataType>&,
        const Element<DataType>&,
        int, int,
        Quadrature<DataType>&,
        Quadrature<DataType>&
        > IFQuadratureSelectionFun;

        DGGlobalAssembler ( );

        ///
        /// \brief Assemble global matrix for a VectorSpace.
        ///
        /// \details Assembles a matrix
        /// \f$A_{ij} = \sum_{K \in \mathcal T_h}\int_{\partial K}{f(x,\varphi_i, \varphi_j)dx}\f$
        /// over the interfaces defined by the mesh
        /// associated to a VectorSpace. The integrand is defined through
        /// a InterfaceMatrixAssemblyFun, which should return the locally
        /// assembled matrix for each element.
        ///
        /// \param[in] space           the VectorSpace for which the assembly is performed
        /// \param[in] local_asm       function or functor that performs local matrix assembly
        /// \param[out] global_matrix  the assembled matrix \f$A_{ij}\f$
        /// \see concept_assembly
        ///
        void assemble_interface_matrix ( const VectorSpace<DataType>& space,
                                         InterfaceMatrixAssemblyFun local_asm,
                                         typename GlobalAssembler<DataType>::GlobalMatrix& matrix ) const;

        ///
        /// \brief Assemble global vector for a VectorSpace.
        ///
        /// \details Assembles a vector
        /// \f$b_i = \sum_{K \in \mathcal T_h} \int_{\partial K}{f(x, \varphi_i)dx}\f$
        /// over the interfaces defined by the mesh
        /// associated to a VectorSpace. The integrand is defined through
        /// a InterfaceVectorAssemblyFun, which should return the locally
        /// assembled vector for each element.
        ///
        /// \param[in] space      the VectorSpace for which the assembly is performed
        /// \param[in] local_asm  function or functor that performs local vector assembly
        /// \param[out] global_vector  the assembled vector \f$b_i\f$
        /// \see concept_assembly
        ///
        void assemble_interface_vector ( const VectorSpace<DataType>& space,
                                         InterfaceVectorAssemblyFun local_asm,
                                         typename GlobalAssembler<DataType>::GlobalVector& vec ) const;

        ///
        /// \brief Assemble global value for a VectorSpace.
        ///
        /// \details Assembles a scalar
        /// \f$v_K = \int_{\partial K}{f(x, \varphi_i)dx}, \ K \in \mathcal T_h\f$
        /// over the interfaces defined by the mesh
        /// associated to a VectorSpace. The integrand is defined through
        /// a InterfaceScalarAssemblyFun, which should return the locally
        /// assembled value for each element.
        ///
        /// \param[in] space      the VectorSpace for which the assembly is performed
        /// \param[in] local_asm  function or functor that performs local vector assembly
        /// \param[out] values  the assembled values \f$v_K\f$
        /// \see concept_assembly
        ///
        void assemble_interface_scalar ( const VectorSpace<DataType>& space,
                                         InterfaceScalarAssemblyFun local_asm,
                                         std::vector<DataType>& values ) const;

        ///
        /// \brief Assemble global value for a VectorSpace.
        ///
        /// \details Assembles a scalar
        /// \f$v_K = \int_{\partial K}{f(x, \varphi_i)dx}, \ K \in \mathcal T_h\f$
        /// over the interfaces defined by the mesh
        /// associated to a VectorSpace. The integrand is defined through
        /// a InterfaceScalarAssemblyFun, which should return the locally
        /// assembled value for each interface.
        ///
        /// \param[in] space      the VectorSpace for which the assembly is performed
        /// \param[in] local_asm  function or functor that performs local vector assembly
        /// \param[out] values  the assembled values \f$v_K\f$ correctly distributed to cells
        /// \see concept_assembly
        ///
        void assemble_interface_scalar_cells ( const VectorSpace<DataType>& space,
                                               InterfaceScalarAssemblyFun local_asm,
                                               std::vector<DataType>& values ) const;

        ///
        /// \brief Assemble global value for a VectorSpace.
        ///
        /// \details Assembles multiple scalars
        /// \f$v_K = \int_{\partial K}{f(x, \varphi_i)dx}, \ K \in \mathcal T_h\f$
        /// over the interfaces defined by the mesh
        /// associated to a VectorSpace. The integrand is defined through
        /// a InterfaceMultipleScalarAssemblyFun, which should return the locally
        /// assembled values for each interface.
        ///
        /// \param[in] space      the VectorSpace for which the assembly is performed
        /// \param[in] local_asm  function or functor that performs local vector assembly
        /// \param[in] num_scalars dimension of integrand
        /// \param[out] values  the assembled values \f$v_K\f$ correctly distributed to cells
        /// \see concept_assembly
        ///
        void assemble_interface_multiple_scalar_cells ( const VectorSpace<DataType>& space,
                                                        InterfaceMultipleScalarAssemblyFun local_asm,
                                                        const int num_scalars,
                                                        std::vector<typename GlobalAssembler<DataType>::LocalVector>& values ) const;

        ///
        /// \brief Add contributions of an interface scalar assembly to a vector of cell values in
        /// a naive way. 
        /// This functionality is, e.g., needed in error estimators
        void distribute_interface_to_cell_values_naive ( const VectorSpace<DataType>& space,
                                                         std::vector<DataType> &cell_values,
                                                         const std::vector<DataType> &interface_values ) const;

        ///
        /// \brief Set the function used to determine which quadrature rule should be used.
        /// for interface quadrature.
        ///
        /// \details The choice of quadrature rules can be controlled
        /// by providing a IFQuadratureSelectionFunction. This function
        /// or functor will be called before the local assembly is
        /// performed on each element. By default, the
        /// DefaultInterfaceQuadratureSelection function is used.
        ///
        /// \param[in] q_select   new quadrature selection function
        void set_interface_quadrature_selection_fun ( IFQuadratureSelectionFun q_select );

      private:
        // special construction of sparsity structure, since dofs do not
        // automatically couple over cell interfaces.
        virtual void compute_sparsity_structure_impl ( const VectorSpace<DataType>& space,
                                                       SparsityStructure& sparsity,
                                                       std::vector<std::vector<bool> > *coupling_vars ) const;

        IFQuadratureSelectionFun if_q_select_;
    };

}

#endif
