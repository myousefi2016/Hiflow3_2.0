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

#include "assembly.h"

#include "generic_assembly_algorithm.h"

#include "assembly/integration_utilities.h"
#include "common/pointers.h"
#include "common/sort_permutation.h"
#include "mesh/attributes.h"
#include "mesh/types.h"
#include "space/element.h"
#include "space/vector_space.h"

namespace hiflow
{

    //////////////// StandardAssembly helper functions ////////////////

    template<template<class, class> class AlgorithmType, class DataType>
    class StandardScalarAssembly
    : public AssemblyAlgorithmBase<AlgorithmType,
    StandardScalarAssembly<AlgorithmType, DataType>, DataType >
    {
      public:
        typedef DataType LocalObjectType;
        typedef hiflow::Quadrature<DataType> QuadratureType;
        typedef AssemblyAlgorithmBase<AlgorithmType, StandardScalarAssembly<AlgorithmType, DataType>, DataType > Base;

        // Name resolution does not manage to get the base members, so
        // we must do it ourselves.
        using Base::space_;
        using Base::curr_;
        using Base::traversal_;
        using Base::elem_;
        using Base::has_next;

        StandardScalarAssembly ( const VectorSpace<DataType>& space, std::vector<DataType>& values )
        : Base ( space ),
        values_ ( values )
        {
            const int num_elements = this->traversal_.size ( );
            this->values_.resize ( num_elements, 0. );

            this->remove_non_local_elements ( );
            sort_elements ( space, this->traversal_ );
        }

        /// The has_next() and next() functions are overloaded in
        /// order to skip elements that do not belong to the local
        /// subdomain. This is done by setting the corresponding
        /// entries in the traversal_ array to -1, and later skipping
        /// those items.

        const Element<DataType>& next ( )
        {
            assert ( this->has_next ( ) );

            this->elem_ = Element<DataType>( this->space_, this->traversal_[this->curr_] );

            ++( this->curr_ );

            return this->elem_;
        }

        void add ( const Element<DataType>& element, const LocalObjectType& local_val )
        {
            this->values_[element.get_cell_index ( )] += local_val;
        }

      private:
        /// Remove non_local elements

        void remove_non_local_elements ( )
        {

            const mesh::Mesh& mesh = this->space_.mesh ( );

            if ( !mesh.has_attribute ( "_remote_index_", mesh.tdim ( ) ) )
            {
                // If the "_remote_index_" attribute does not exist, we
                // assume that there are no ghost cells.
                return;
            }

            int remote_index;
            int index;
            std::vector<int> traversal_tmp;
            traversal_tmp.reserve ( mesh.num_entities ( mesh.tdim ( ) ) );

            for ( mesh::EntityIterator it_cell = mesh.begin ( mesh.tdim ( ) ), e_it_cell = mesh.end ( mesh.tdim ( ) );
                  it_cell != e_it_cell;
                  ++it_cell )
            {
                // test if cell on subdomain
                it_cell->get ( "_remote_index_", &remote_index );
                if ( remote_index == -1 )
                {
                    index = it_cell->index ( );
                    traversal_tmp.push_back ( index );
                }
            }
            this->traversal_ = traversal_tmp;
        }

        std::vector<DataType>& values_;
    };

    template<template<class, class> class AlgorithmType, class DataType>
    class StandardMultipleScalarAssembly
    : public AssemblyAlgorithmBase<AlgorithmType,
    StandardMultipleScalarAssembly<AlgorithmType, DataType>, DataType >
    {
      public:
        typedef std::vector<DataType> LocalObjectType;
        typedef hiflow::Quadrature<DataType> QuadratureType;
        typedef AssemblyAlgorithmBase<AlgorithmType, StandardMultipleScalarAssembly<AlgorithmType, DataType>, DataType > Base;

        // Name resolution does not manage to get the base members, so
        // we must do it ourselves.
        using Base::space_;
        using Base::curr_;
        using Base::traversal_;
        using Base::elem_;
        using Base::has_next;

        StandardMultipleScalarAssembly ( const VectorSpace<DataType>& space, std::vector< std::vector<DataType> >& values, const int num_scalars )
        : Base ( space ),
        values_ ( values ),
        num_scalars_ ( num_scalars )
        {
            const int num_elements = this->traversal_.size ( );
            this->values_.resize ( num_elements );
            for ( int l = 0; l < num_elements; ++l )
            {
                this->values_[l].resize ( this->num_scalars_, 0. );
            }
            this->remove_non_local_elements ( );
            sort_elements ( space, this->traversal_ );
        }

        /// The has_next() and next() functions are overloaded in
        /// order to skip elements that do not belong to the local
        /// subdomain. This is done by setting the corresponding
        /// entries in the traversal_ array to -1, and later skipping
        /// those items.

        const Element<DataType>& next ( )
        {
            assert ( this->has_next ( ) );

            this->elem_ = Element<DataType>( this->space_, this->traversal_[this->curr_] );

            ++( this->curr_ );

            return this->elem_;
        }

        void add ( const Element<DataType>& element, const LocalObjectType& local_val )
        {
            for ( int l = 0; l < this->num_scalars_; ++l )
            {
                this->values_[element.get_cell_index ( )][l] += local_val[l];
            }
        }

      private:
        /// Remove non_local elements

        void remove_non_local_elements ( )
        {

            const mesh::Mesh& mesh = this->space_.mesh ( );

            if ( !mesh.has_attribute ( "_remote_index_", mesh.tdim ( ) ) )
            {
                // If the "_remote_index_" attribute does not exist, we
                // assume that there are no ghost cells.
                return;
            }

            int remote_index;
            int index;
            std::vector<int> traversal_tmp;
            traversal_tmp.reserve ( mesh.num_entities ( mesh.tdim ( ) ) );

            for ( mesh::EntityIterator it_cell = mesh.begin ( mesh.tdim ( ) ), e_it_cell = mesh.end ( mesh.tdim ( ) );
                  it_cell != e_it_cell;
                  ++it_cell )
            {
                // test if cell on subdomain
                it_cell->get ( "_remote_index_", &remote_index );
                if ( remote_index == -1 )
                {
                    index = it_cell->index ( );
                    traversal_tmp.push_back ( index );
                }
            }
            this->traversal_ = traversal_tmp;
        }
        int num_scalars_;
        std::vector< std::vector<DataType> >& values_;
    };

    template<template<class, class> class AlgorithmType, class DataType>
    class StandardVectorAssembly
    : public AssemblyAlgorithmBase< AlgorithmType, StandardVectorAssembly<AlgorithmType, DataType>, DataType >
    {
      public:
        typedef std::vector<DataType> LocalObjectType;
        typedef hiflow::Quadrature<DataType> QuadratureType;
        typedef AssemblyAlgorithmBase< AlgorithmType, StandardVectorAssembly<AlgorithmType, DataType>, DataType > Base;

        // Name resolution does not manage to get the base members, so
        // we must do it ourselves.
        using Base::space_;
        using Base::dof_;
        using Base::traversal_;

        StandardVectorAssembly ( const VectorSpace<DataType>& space, typename GlobalAssembler<DataType>::GlobalVector& vec )
        : Base ( space ), vector_ ( vec )
        {

            sort_elements ( space, this->traversal_ );
        }

        void add ( const Element<DataType>& element, const LocalObjectType& local_vec )
        {
            const int num_dofs = this->dof_.size ( );

            std::vector<int> dofs_sort_permutation;

            // get permutation for sorting dofs
            sortingPermutation ( this->dof_, dofs_sort_permutation );

            // create row array
            std::vector<int> row_indices;
            row_indices.reserve ( num_dofs );

            LocalObjectType local_vec_sorted;
            local_vec_sorted.reserve ( num_dofs );

            for ( size_t i = 0; i != num_dofs; ++i )
            {
                const int dof_sort_perm = dofs_sort_permutation[i];
                const int dof_ind = this->dof_[dof_sort_perm];
                if ( this->space_.dof ( ).is_dof_on_sd ( dof_ind ) )
                {
                    row_indices.push_back ( dof_ind );
                    local_vec_sorted.push_back ( local_vec[dof_sort_perm] );
                }
            }

            // Add local to global vector
            if ( row_indices.size ( ) > 0 )
            {
                this->vector_.Add ( vec2ptr ( row_indices ), row_indices.size ( ), vec2ptr ( local_vec_sorted ) );
            }
        }

      private:
        typename GlobalAssembler<DataType>::GlobalVector& vector_;

    };

    template<template<class, class> class AlgorithmType, class DataType>
    class StandardMatrixAssembly
    : public AssemblyAlgorithmBase<AlgorithmType, StandardMatrixAssembly<AlgorithmType, DataType>, DataType >
    {
      public:
        typedef la::SeqDenseMatrix<DataType> LocalObjectType;
        typedef Quadrature<DataType> QuadratureType;
        typedef AssemblyAlgorithmBase<AlgorithmType, StandardMatrixAssembly<AlgorithmType, DataType>, DataType > Base;

        // Name resolution does not manage to get the base members, so
        // we must do it ourselves.
        using Base::space_;
        using Base::dof_;
        using Base::traversal_;

        StandardMatrixAssembly ( const VectorSpace<DataType>& space, typename GlobalAssembler<DataType>::GlobalMatrix& matrix )
        : Base ( space ), matrix_ ( matrix )
        {
            sort_elements ( space, this->traversal_ );
        }

        void add ( const Element<DataType>& element, const LocalObjectType& local_mat )
        {
            const int num_dofs = this->dof_.size ( );

            std::vector<int> dofs_sort_permutation;
            std::vector<int> dofs_sorted ( num_dofs );

            // get permutation for sorting dofs
            sortingPermutation ( this->dof_, dofs_sort_permutation );

            // fill sorted dof array
            for ( size_t i = 0; i != num_dofs; ++i )
            {
                dofs_sorted[i] = this->dof_[dofs_sort_permutation[i]];
            }

            // create row array
            std::vector<int> row_indices;
            std::vector<int> row_permutation;
            row_indices.reserve ( num_dofs );
            row_permutation.reserve ( num_dofs );
            for ( size_t i = 0; i != num_dofs; ++i )
            {
                const int dof_sort_perm = dofs_sort_permutation[i];
                const int dof_ind = this->dof_[dof_sort_perm];
                if ( this->space_.dof ( ).is_dof_on_sd ( dof_ind ) )
                {
                    row_indices.push_back ( dof_ind );
                    row_permutation.push_back ( dof_sort_perm );
                }
            }

            // fill reduced and sorted local matrix
            LocalObjectType local_mat_sorted_reduced;
            if ( row_indices.size ( ) > 0 && num_dofs > 0 )
            {
                local_mat_sorted_reduced.Resize ( row_indices.size ( ), num_dofs );
                for ( size_t i = 0, i_e = row_indices.size ( ); i != i_e; ++i )
                {
                    const int row_ind = row_permutation[i];
                    for ( size_t j = 0, j_e = num_dofs; j != j_e; ++j )
                    {
                        local_mat_sorted_reduced ( i, j ) = local_mat ( row_ind, dofs_sort_permutation[j] );
                    }
                }

                // Add local to global matrix
                this->matrix_.Add ( vec2ptr ( row_indices ), row_indices.size ( ), vec2ptr ( dofs_sorted ), dofs_sorted.size ( ), &local_mat_sorted_reduced ( 0, 0 ) );
            }
        }

      private:
        typename GlobalAssembler<DataType>::GlobalMatrix& matrix_;

    };

    template<template<class, class> class AlgorithmType, class DataType>
    class StandardBoundaryScalarAssembly
    : public AssemblyAlgorithmBase<AlgorithmType,
    StandardBoundaryScalarAssembly<AlgorithmType, DataType>, DataType >
    {
      public:
        typedef std::vector<DataType> LocalObjectType;
        typedef hiflow::Quadrature<DataType> QuadratureType;
        typedef AssemblyAlgorithmBase<AlgorithmType, StandardBoundaryScalarAssembly<AlgorithmType, DataType>, DataType > Base;

        // Name resolution does not manage to get the base members, so
        // we must do it ourselves.
        using Base::space_;
        using Base::curr_;
        using Base::traversal_;
        using Base::elem_;
        using Base::has_next;

        StandardBoundaryScalarAssembly ( const VectorSpace<DataType>& space, std::vector<DataType>& values )
        : Base ( space ),
        values_ ( values )
        {
            remove_non_local_elements ( );
            sort_elements ( space, this->traversal_ );
        }

        /// The has_next() and next() functions are overloaded in
        /// order to skip elements that do not belong to the local
        /// subdomain. This is done by setting the corresponding
        /// entries in the traversal_ array to -1, and later skipping
        /// those items.

        const Element<DataType>& next ( )
        {
            assert ( this->has_next ( ) );

            this->elem_ = Element<DataType>( this->space_, this->traversal_[this->curr_] );

            ++( this->curr_ );

            return this->elem_;
        }

        void add ( const Element<DataType>& element, const LocalObjectType& local_val )
        {
            mesh::TDim tdim = this->space_.mesh ( ).tdim ( );
            mesh::IncidentEntityIterator iter = element.get_cell ( ).begin_incident ( tdim - 1 );
            mesh::IncidentEntityIterator end = element.get_cell ( ).end_incident ( tdim - 1 );
            int facet_number = 0;
            for (; iter != end; iter++ )
            {

                this->values_[iter->id ( )] += local_val[facet_number];
                ++facet_number;
            }
        }

        void reset ( typename GlobalAssembler<DataType>::LocalVector& local_vec )
        {
            mesh::TDim tdim = this->space_.mesh ( ).tdim ( );
            local_vec.clear ( );
            local_vec.resize ( this->elem_.get_cell ( ).num_incident_entities ( tdim - 1 ), 0. );
        }

      private:
        /// Remove non_local elements

        void remove_non_local_elements ( )
        {

            const mesh::Mesh& mesh = this->space_.mesh ( );

            if ( !mesh.has_attribute ( "_remote_index_", mesh.tdim ( ) ) )
            {
                // If the "_remote_index_" attribute does not exist, we
                // assume that there are no ghost cells.
                return;
            }

            int remote_index;
            int index;
            std::vector<int> traversal_tmp;
            traversal_tmp.reserve ( mesh.num_entities ( mesh.tdim ( ) ) );

            for ( mesh::EntityIterator it_cell = mesh.begin ( mesh.tdim ( ) ), e_it_cell = mesh.end ( mesh.tdim ( ) );
                  it_cell != e_it_cell;
                  ++it_cell )
            {
                // test if cell on subdomain
                it_cell->get ( "_remote_index_", &remote_index );
                if ( remote_index == -1 )
                {
                    index = it_cell->index ( );
                    traversal_tmp.push_back ( index );
                }
            }
            this->traversal_ = traversal_tmp;
        }

        std::vector<DataType>& values_;
    };
    //////////////// end helper functions ////////////////

    //////////////// Implementation of StandardGlobalAssembler ////////////////

    template<class DataType>
    void StandardGlobalAssembler<DataType>::compute_sparsity_structure_impl (
                                                                              const VectorSpace<DataType>& space,
                                                                              SparsityStructure& sparsity,
                                                                              std::vector<std::vector<bool> > *coupling_vars ) const
    {
        // Assert correct size of coupling_vars

        assert ( coupling_vars->size ( ) == space.get_nb_var ( ) );
        for ( size_t i = 0; i < space.get_nb_var ( ); ++i )
        {
            assert ( ( *coupling_vars )[i].size ( ) == space.get_nb_var ( ) );
        }

        // Call function from integration_utilities.h
        // TODO (staffan) Refactor this in future.
        InitStructure ( space,
                        &sparsity.diagonal_rows,
                        &sparsity.diagonal_cols,
                        &sparsity.off_diagonal_rows,
                        &sparsity.off_diagonal_cols,
                        coupling_vars );
    }

    template<class DataType>
    void StandardGlobalAssembler<DataType>::assemble_scalar_impl (
                                                                   const VectorSpace<DataType>& space,
                                                                   typename GlobalAssembler<DataType>::ScalarAssemblyFunction local_asm,
                                                                   std::vector<DataType>& vec,
                                                                   typename GlobalAssembler<DataType>::QuadratureSelectionFunction q_select ) const
    {
        StandardScalarAssembly<InteriorAssemblyAlgorithm, DataType> assembly ( space, vec );
        assembly.assemble ( local_asm, q_select );
    }

    template<class DataType>
    void StandardGlobalAssembler<DataType>::assemble_multiple_scalar_impl (
                                                                            const VectorSpace<DataType>& space,
                                                                            typename GlobalAssembler<DataType>::MultipleScalarAssemblyFunction local_asm,
                                                                            const int num_scalars,
                                                                            std::vector< std::vector<DataType> >& vec,
                                                                            typename GlobalAssembler<DataType>::QuadratureSelectionFunction q_select ) const
    {
        StandardMultipleScalarAssembly<InteriorAssemblyAlgorithm, DataType> assembly ( space, vec, num_scalars );
        assembly.assemble ( local_asm, q_select );
    }

    template<class DataType>
    void StandardGlobalAssembler<DataType>::assemble_vector_impl (
                                                                   const VectorSpace<DataType>& space,
                                                                   typename GlobalAssembler<DataType>::VectorAssemblyFunction local_asm,
                                                                   typename GlobalAssembler<DataType>::GlobalVector& vec,
                                                                   typename GlobalAssembler<DataType>::QuadratureSelectionFunction q_select ) const
    {
        StandardVectorAssembly<InteriorAssemblyAlgorithm, DataType> assembly ( space, vec );
        assembly.assemble ( local_asm, q_select );
    }

    template<class DataType>
    void StandardGlobalAssembler<DataType>::assemble_matrix_impl (
                                                                   const VectorSpace<DataType>& space,
                                                                   typename GlobalAssembler<DataType>::MatrixAssemblyFunction local_asm,
                                                                   typename GlobalAssembler<DataType>::GlobalMatrix& mat,
                                                                   typename GlobalAssembler<DataType>::QuadratureSelectionFunction q_select ) const
    {
        StandardMatrixAssembly<InteriorAssemblyAlgorithm, DataType> assembly ( space, mat );
        assembly.assemble ( local_asm, q_select );
    }

    template<class DataType>
    void StandardGlobalAssembler<DataType>::assemble_scalar_boundary_impl (
                                                                            const VectorSpace<DataType>& space,
                                                                            typename GlobalAssembler<DataType>::BoundaryScalarAssemblyFunction local_asm,
                                                                            std::vector<DataType>& vec,
                                                                            typename GlobalAssembler<DataType>::FacetQuadratureSelectionFunction fq_select ) const
    {
        // TODO: how should vec be defined -> all facets or only boundary facets?
        // If the latter, then what ordering should be used?
        StandardBoundaryScalarAssembly<BoundaryAssemblyAlgorithm, DataType> assembly ( space, vec );
        assembly.assemble ( local_asm, fq_select );
    }

    template<class DataType>
    void StandardGlobalAssembler<DataType>::assemble_vector_boundary_impl (
                                                                            const VectorSpace<DataType>& space,
                                                                            typename GlobalAssembler<DataType>::BoundaryVectorAssemblyFunction local_asm,
                                                                            typename GlobalAssembler<DataType>::GlobalVector& vec,
                                                                            typename GlobalAssembler<DataType>::FacetQuadratureSelectionFunction fq_select ) const
    {
        StandardVectorAssembly<BoundaryAssemblyAlgorithm, DataType> assembly ( space, vec );
        assembly.assemble ( local_asm, fq_select );
    }

    template<class DataType>
    void StandardGlobalAssembler<DataType>::assemble_matrix_boundary_impl (
                                                                            const VectorSpace<DataType>& space,
                                                                            typename GlobalAssembler<DataType>::BoundaryMatrixAssemblyFunction local_asm,
                                                                            typename GlobalAssembler<DataType>::GlobalMatrix& mat,
                                                                            typename GlobalAssembler<DataType>::FacetQuadratureSelectionFunction fq_select ) const
    {
        StandardMatrixAssembly<BoundaryAssemblyAlgorithm, DataType> assembly ( space, mat );
        assembly.assemble ( local_asm, fq_select );
    }

    template class StandardGlobalAssembler<double>;
    template class StandardGlobalAssembler<float>;
}
