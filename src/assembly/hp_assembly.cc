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
#include "common/log.h"
#include "common/pointers.h"
#include "common/sorted_array.h"
#include "dof/dof_interpolation.h"
#include "mesh/types.h"
#include "space/element.h"
#include "space/vector_space.h"
#include "common/log.h"
#include <utility>
#include <vector>

const int DEBUG_LEVEL = 0;
//#define OCTAVE_OUTPUT

namespace hiflow
{
    typedef std::vector< std::pair<int, double> >::const_iterator ConstraintIterator;

    //using doffem::DofInterpolation;
    using namespace doffem;

    template<class DataType>
    void compute_hp_sparsity_structure ( const VectorSpace<DataType>& space, SparsityStructure& sparsity, std::vector<std::vector<bool> > *coupling_vars )
    {

        // Assert correct size of coupling_vars

        assert ( coupling_vars->size ( ) == space.get_nb_var ( ) );
        for ( size_t i = 0, e_i = space.get_nb_var ( ); i != e_i; ++i )
        {
            assert ( ( *coupling_vars )[i].size ( ) == space.get_nb_var ( ) );
        }

        // TODO: refactor function to avoid all the repetitions

        typedef typename VectorSpace<DataType>::MeshEntityIterator CellIterator;
        typedef std::vector<int>::const_iterator DofIterator;
        const mesh::Mesh& mesh = space.mesh ( );
        const mesh::TDim tdim = mesh.tdim ( );

        std::vector<int> dofs_test, dofs_trial;

        const DofInterpolation& interpolation = space.dof ( ).dof_interpolation ( );
        const int num_total_dofs = space.dof ( ).ndofs_on_sd ( space.dof ( ).get_my_subdomain ( ) );
        int local_row_dof;

        // NB: We assume that unconstrained dofs are numbered before
        // constrained dofs, in order to be able to use a vector here.
        std::vector< SortedArray<int> > diagonal_couplings ( num_total_dofs );
        std::vector< SortedArray<int> > off_diagonal_couplings ( num_total_dofs );

        // Loop over all cells
        for ( CellIterator cell_it = mesh.begin ( tdim ), cell_end = mesh.end ( tdim );
              cell_it != cell_end; ++cell_it )
        {

            // loop over test variables
            for ( size_t test_var = 0, e_test_var = space.get_nb_var ( ); test_var != e_test_var; ++test_var )
            {
                // Get dof id:s on cell
                space.GetDofIndices ( test_var, *cell_it, &dofs_test );

                // Loop over rows corresponding to local dofs.
                for ( DofIterator it_i = dofs_test.begin ( ), end_i = dofs_test.end ( ); it_i != end_i; ++it_i )
                {

                    // search current dof it_i in DofInterpolation map 
                    DofInterpolation::const_iterator dof_i = interpolation.find ( *it_i );

                    if ( dof_i != interpolation.end ( ) )
                    {
                        // Case A: dofs_test[*it_i] (row) constrained

                        // Loop over all interpolating dofs of current dof it_i
                        for ( ConstraintIterator ci_it = dof_i->second.begin ( ), ci_end = dof_i->second.end ( ); ci_it != ci_end; ++ci_it )
                        {

                            // skip rows that are not on our sub-domain
                            if ( space.dof ( ).is_dof_on_sd ( ci_it->first ) )
                            {

                                // get row index to use for insertion
                                space.dof ( ).global2local ( ci_it->first, &local_row_dof );

                                // loop over trial variables
                                for ( size_t trial_var = 0, e_trial_var = space.get_nb_var ( ); trial_var != e_trial_var; ++trial_var )
                                {

                                    // check if coupling exists
                                    if ( ( *coupling_vars )[test_var][trial_var] )
                                    {

                                        // Get dof id:s on cell
                                        space.GetDofIndices ( trial_var, *cell_it, &dofs_trial );

                                        // Loop over columns corresponding to local dofs.
                                        for ( DofIterator it_j = dofs_trial.begin ( ), end_j = dofs_trial.end ( );
                                              it_j != end_j; ++it_j )
                                        {

                                            // search current dof it_j in DofInterpolation map 
                                            DofInterpolation::const_iterator dof_j = interpolation.find ( *it_j );

                                            if ( dof_j != interpolation.end ( ) )
                                            {
                                                // Case A1: dofs_trial[*it_j] (column) constrained

                                                // Loop over all interpolating dofs of current dof it_j
                                                for ( ConstraintIterator cj_it = dof_j->second.begin ( ), cj_end = dof_j->second.end ( ); cj_it != cj_end; ++cj_it )
                                                {
                                                    // determine target for insertion                                                    
                                                    if ( space.dof ( ).is_dof_on_sd ( cj_it->first ) )
                                                    {
                                                        diagonal_couplings[local_row_dof].find_insert ( cj_it->first );
                                                    }
                                                    else
                                                    {
                                                        off_diagonal_couplings[local_row_dof].find_insert ( cj_it->first );
                                                    }

                                                    LOG_DEBUG ( 2, "[" << space.dof ( ).get_my_subdomain ( ) << "]   Constrained row = " << local_row_dof
                                                                << ", constrained col = " << cj_it->first
                                                                << ", diagonal ? " << space.dof ( ).is_dof_on_sd ( cj_it->first ) );

                                                }
                                            }
                                            else
                                            {
                                                // Case A2: dofs_trial[*it_j] (column) unconstrained
                                                // determine target for insertion
                                                if ( space.dof ( ).is_dof_on_sd ( *it_j ) )
                                                {
                                                    diagonal_couplings[local_row_dof].find_insert ( *it_j );
                                                }
                                                else
                                                {
                                                    off_diagonal_couplings[local_row_dof].find_insert ( *it_j );
                                                }

                                                LOG_DEBUG ( 2, "[" << space.dof ( ).get_my_subdomain ( ) << "]   Constrained row = " << local_row_dof
                                                            << ", unconstrained col = " << *it_j
                                                            << ", diagonal ? " << space.dof ( ).is_dof_on_sd ( *it_j ) );
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    else
                    {
                        // Case B: dofs_test[*it_i] (row) unconstrained

                        // skip rows that are not on our sub-domain
                        if ( space.dof ( ).is_dof_on_sd ( *it_i ) )
                        {

                            space.dof ( ).global2local ( *it_i, &local_row_dof );

                            // loop over trial variables
                            for ( size_t trial_var = 0, e_trial_var = space.get_nb_var ( ); trial_var != e_trial_var; ++trial_var )
                            {

                                // check if coupling exists
                                if ( ( *coupling_vars )[test_var][trial_var] )
                                {

                                    // Get dof id:s on cell
                                    space.GetDofIndices ( trial_var, *cell_it, &dofs_trial );

                                    // Loop over columns corresponding to local dofs.
                                    for ( DofIterator it_j = dofs_trial.begin ( ), end_j = dofs_trial.end ( ); it_j != end_j; ++it_j )
                                    {

                                        // search current dof it_j in DofInterpolation map
                                        DofInterpolation::const_iterator dof_j = interpolation.find ( *it_j );

                                        if ( dof_j != interpolation.end ( ) )
                                        {
                                            // Case B1: dofs_trial[*it_j] (column) constrained

                                            // Loop over all interpolating dofs of current dof it_j
                                            for ( ConstraintIterator cj_it = dof_j->second.begin ( ), cj_end = dof_j->second.end ( ); cj_it != cj_end; ++cj_it )
                                            {
                                                // determine target for insertion
                                                // -> diagonal or off-diagonal
                                                if ( space.dof ( ).is_dof_on_sd ( cj_it->first ) )
                                                {
                                                    diagonal_couplings[local_row_dof].find_insert ( cj_it->first );
                                                }
                                                else
                                                {
                                                    off_diagonal_couplings[local_row_dof].find_insert ( cj_it->first );
                                                }

                                                LOG_DEBUG ( 2, "[" << space.dof ( ).get_my_subdomain ( ) << "] Unconstrained row = " << local_row_dof
                                                            << ", constrained col = " << cj_it->first
                                                            << ", diagonal ? " << space.dof ( ).is_dof_on_sd ( cj_it->first ) );
                                            }
                                        }
                                        else
                                        {
                                            // Case B2: dofs_trial[*it_j] (column) unconstrained
                                            // determine target for insertion
                                            // -> diagonal or off-diagonal
                                            if ( space.dof ( ).is_dof_on_sd ( *it_j ) )
                                            {
                                                diagonal_couplings[local_row_dof].find_insert ( *it_j );
                                            }
                                            else
                                            {
                                                off_diagonal_couplings[local_row_dof].find_insert ( *it_j );
                                            }

                                            LOG_DEBUG ( 3, "[" << space.dof ( ).get_my_subdomain ( ) << "] Unconstrained row = " << local_row_dof
                                                        << ", unconstrained col = " << *it_j
                                                        << ", diagonal ? " << space.dof ( ).is_dof_on_sd ( *it_j ) );
                                        }
                                    }
                                }
                            }
                        }
                    }
                    // Add diagonal entry
                    if ( space.dof ( ).is_dof_on_sd ( *it_i ) )
                    {
                        space.dof ( ).global2local ( *it_i, &local_row_dof );
                        diagonal_couplings[local_row_dof].find_insert ( *it_i );

                        LOG_DEBUG ( 2, "[" << space.dof ( ).get_my_subdomain ( ) << "]   Diagonal row = " << local_row_dof << ", diagonal col = " << local_row_dof );
                    }
                }
            }
        }

        // Compute number of non-zeros for both blocks
        int nnz_diagonal = 0;
        int nnz_off_diagonal = 0;

        for ( size_t i = 0; i != num_total_dofs; ++i )
        {
            nnz_diagonal += diagonal_couplings[i].size ( );
            nnz_off_diagonal += off_diagonal_couplings[i].size ( );
        }

        // Copy into SparsityStructure
        sparsity.diagonal_rows.resize ( nnz_diagonal );
        sparsity.diagonal_cols.resize ( nnz_diagonal );
        sparsity.off_diagonal_rows.resize ( nnz_off_diagonal );
        sparsity.off_diagonal_cols.resize ( nnz_off_diagonal );

        int global_dof_id;
        int index = 0;
        for ( size_t i = 0; i != num_total_dofs; ++i )
        {
            space.dof ( ).local2global ( i, &global_dof_id );

            for ( SortedArray<int>::const_iterator it = diagonal_couplings[i].begin ( ),
                  end = diagonal_couplings[i].end ( ); it != end; ++it )
            {
                sparsity.diagonal_rows[index] = global_dof_id;
                sparsity.diagonal_cols[index] = *it;
                ++index;
            }
        }
        assert ( index == nnz_diagonal );

        index = 0;
        for ( size_t i = 0; i != num_total_dofs; ++i )
        {
            space.dof ( ).local2global ( i, &global_dof_id );

            for ( SortedArray<int>::const_iterator it = off_diagonal_couplings[i].begin ( ),
                  end = off_diagonal_couplings[i].end ( ); it != end; ++it )
            {
                sparsity.off_diagonal_rows[index] = global_dof_id;
                sparsity.off_diagonal_cols[index] = *it;
                ++index;
            }
        }
        assert ( index == nnz_off_diagonal );
    }

    template<class DataType>
    class HpVectorAssembly
    : public AssemblyAlgorithmBase<InteriorAssemblyAlgorithm, HpVectorAssembly<DataType>, DataType>
    {
      public:
        typedef typename GlobalAssembler<DataType>::LocalVector LocalObjectType;
        typedef hiflow::Quadrature<DataType> QuadratureType;

        HpVectorAssembly ( const VectorSpace<DataType>& space, typename GlobalAssembler<DataType>::GlobalVector& vec )
        : AssemblyAlgorithmBase<hiflow::InteriorAssemblyAlgorithm, HpVectorAssembly, DataType>( space ),
        vector_ ( vec ),
        interp_ ( space.dof ( ).dof_interpolation ( ) )
        {
#ifdef OCTAVE_OUTPUT
            octave_.open ( "check_assembly.m", std::ios_base::app );
            octave_.precision ( 16 );
            octave_ << "% Global vector assembly\n";
            octave_ << "b = zeros(" << vector_.size_global ( ) << ", 1);\n";
#endif
            vector_.Zeros ( );
        }

        ~HpVectorAssembly ( )
        {
#ifdef OCTAVE_OUTPUT
            octave_ << "\n\n";
            octave_.close ( );
#endif
        }

        void add ( const Element<DataType>& element, const LocalObjectType& local_vec )
        {
            const int num_dofs = this->dof_.size ( );
            for ( size_t i = 0; i != num_dofs; ++i )
            {
                DofInterpolation::const_iterator it = this->interp_.find ( this->dof_[i] );
                if ( it != this->interp_.end ( ) )
                {
                    // dof[i] is constrained -> add contributions to dependent dofs
                    for ( ConstraintIterator c_it = it->second.begin ( ),
                          c_end = it->second.end ( ); c_it != c_end; ++c_it )
                    {
                        if ( this->space_.dof ( ).is_dof_on_sd ( c_it->first ) )
                        {
                            vector_.Add ( c_it->first, c_it->second * local_vec[i] );
                        }
                    }
                }
                else
                {
                    // dof[i] is unconstrained -> add contribution to this dof
                    if ( this->space_.dof ( ).is_dof_on_sd ( this->dof_[i] ) )
                    {
                        vector_.Add ( this->dof_[i], local_vec[i] );
                    }
                }
            }
#ifdef OCTAVE_OUTPUT
            octave_ << "% Element " << element.get_cell_index ( ) << "\n";
            octave_ << "dof = [" << string_from_range ( dof_.begin ( ), dof_.end ( ) ) << "] + 1;\n"
                    << "b_local = ["
                    << precise_string_from_range ( local_vec.begin ( ), local_vec.end ( ) ) << "]';\n"
                    << "b(dof) += b_local;\n";
#endif
        }

      private:
        typename GlobalAssembler<DataType>::GlobalVector& vector_;
        const DofInterpolation& interp_;

#ifdef OCTAVE_OUTPUT
        std::ofstream octave_;
#endif

    };

    template<class DataType>
    class HpMatrixAssembly
    : public AssemblyAlgorithmBase<InteriorAssemblyAlgorithm, HpMatrixAssembly<DataType>, DataType>
    {
      public:
        typedef la::SeqDenseMatrix<DataType> LocalObjectType;
        typedef Quadrature<DataType> QuadratureType;

        HpMatrixAssembly ( const VectorSpace<DataType>& space, typename GlobalAssembler<DataType>::GlobalMatrix& matrix )
        : AssemblyAlgorithmBase<hiflow::InteriorAssemblyAlgorithm, HpMatrixAssembly, DataType>( space ),
        matrix_ ( matrix ),
        interp_ ( space.dof ( ).dof_interpolation ( ) )
        {
#ifdef OCTAVE_OUTPUT

            octave_.open ( "check_assembly.m", std::ios_base::app );
            octave_.precision ( 16 );
            octave_ << "% ==== Global matrix assembly ====\n";
            octave_ << "A = zeros(" << matrix.nrows_global ( ) << ");\n";
#endif
            matrix_.Zeros ( );
        }

        ~HpMatrixAssembly ( )
        {
            // Set rows of constrained dofs to identity to obtain non-singular
            // matrix
            SortedArray<int> constrained_dofs;
            for ( DofInterpolation::const_iterator it = interp_.begin ( ), end = interp_.end ( );
                  it != end; ++it )
            {

                if ( this->space_.dof ( ).is_dof_on_sd ( it->first ) )
                {
                    constrained_dofs.find_insert ( it->first );
                }
            }

            if ( !constrained_dofs.empty ( ) )
            {
                matrix_.diagonalize_rows ( &constrained_dofs.front ( ), constrained_dofs.size ( ), 1. );
            }
#ifdef OCTAVE_OUTPUT
            // Close octave stream
            octave_ << "\n\n";
            octave_.close ( );
#endif
        }

        void add ( const Element<DataType>& element, const LocalObjectType& local_mat )
        {
            // Assemble into global system.  Only add entries to
            // unconstrained rows and columns, and only if the dof
            // corresponding to the row belongs to the local subdomain.
            const int num_dofs = this->dof_.size ( );
            for ( size_t i = 0; i != num_dofs; ++i )
            {
                DofInterpolation::const_iterator it_i = this->interp_.find ( this->dof_[i] );

                if ( it_i != this->interp_.end ( ) )
                {
                    // dof[i] is constrained -> add contributions to dependent rows
                    for ( ConstraintIterator c_it = it_i->second.begin ( ),
                          c_end = it_i->second.end ( ); c_it != c_end; ++c_it )
                    {
                        if ( this->space_.dof ( ).is_dof_on_sd ( c_it->first ) )
                        {

                            for ( size_t j = 0; j != num_dofs; ++j )
                            {
                                DofInterpolation::const_iterator it_j = this->interp_.find ( this->dof_[j] );

                                if ( it_j != this->interp_.end ( ) )
                                {
                                    // dof[j] is constrained -> add attributions to dependent columns
                                    // TODO: are these not cleared at the end anyway?
                                    for ( ConstraintIterator c2_it = it_j->second.begin ( ),
                                          c2_end = it_j->second.end ( ); c2_it != c2_end; ++c2_it )
                                    {
                                        matrix_.Add ( c_it->first,
                                                      c2_it->first,
                                                      c_it->second * c2_it->second * local_mat ( i, j ) );
                                    }
                                }
                                else
                                {
                                    // dof[j] unconstrained -> add contribution to dof[j] column
                                    matrix_.Add ( c_it->first,
                                                  this->dof_[j],
                                                  c_it->second * local_mat ( i, j ) );
                                }
                            }
                        }
                    }
                }
                else
                {
                    // dof[i] is unconstrained
                    if ( this->space_.dof ( ).is_dof_on_sd ( this->dof_[i] ) )
                    {
                        for ( size_t j = 0; j != num_dofs; ++j )
                        {
                            DofInterpolation::const_iterator it_j = this->interp_.find ( this->dof_[j] );
                            if ( it_j != this->interp_.end ( ) )
                            {
                                for ( ConstraintIterator c_it = it_j->second.begin ( ),
                                      c_end = it_j->second.end ( ); c_it != c_end; ++c_it )
                                {
                                    // dof[j] is constrained -> add attributions to dependent columns
                                    matrix_.Add ( this->dof_[i],
                                                  c_it->first,
                                                  c_it->second * local_mat ( i, j ) );
                                }
                            }
                            else
                            {
                                // dof[j] unconstrained - assemble normally
                                matrix_.Add ( this->dof_[i], this->dof_[j], local_mat ( i, j ) );
                            }
                        }
                    }
                }
            }
#ifdef OCTAVE_OUTPUT
            octave_ << "% Element " << element.get_cell_index ( ) << "\n";
            octave_ << "dof = [" << string_from_range ( dof_.begin ( ), dof_.end ( ) ) << "] + 1;\n"
                    << "A_local = " << local_mat << ";\n"
                    << "A(dof, dof) += A_local;\n";
#endif
        }

      private:
        typename GlobalAssembler<DataType>::GlobalMatrix& matrix_;
        const DofInterpolation& interp_;
#ifdef OCTAVE_OUTPUT
        std::ofstream octave_;
#endif
    };

    //////////////////////////////////////////////////////////////////////////////////
    //////////////// Implementation of HpFemAssembler ////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    template<class DataType>
    void HpFemAssembler<DataType>::compute_sparsity_structure_impl ( const VectorSpace<DataType>& space,
                                                                     SparsityStructure& sparsity,
                                                                     std::vector<std::vector<bool> > *coupling_vars ) const
    {
        // Assert correct size of coupling_vars

        assert ( coupling_vars->size ( ) == space.get_nb_var ( ) );
        for ( size_t i = 0, e_i = space.get_nb_var ( ); i != e_i; ++i )
        {
            assert ( ( *coupling_vars )[i].size ( ) == space.get_nb_var ( ) );
        }

        compute_hp_sparsity_structure ( space, sparsity, coupling_vars );
    }

    template<class DataType>
    void HpFemAssembler<DataType>::assemble_scalar_impl ( const VectorSpace<DataType>& space,
                                                          typename GlobalAssembler<DataType>::ScalarAssemblyFunction local_asm,
                                                          std::vector<DataType>& vec,
                                                          typename GlobalAssembler<DataType>::QuadratureSelectionFunction q_select ) const
    {
        // Exactly same implementation as Standard Assembly, so we delegate.
        // TODO(staffan): maybe we should derive from StandardGlobalAssembler instead?
        StandardGlobalAssembler<DataType> assembly;
        assembly.set_quadrature_selection_function ( q_select );
        assembly.assemble_scalar ( space, local_asm, vec );
    }

    template<class DataType>
    void HpFemAssembler<DataType>::assemble_multiple_scalar_impl ( const VectorSpace<DataType>& space,
                                                                   typename GlobalAssembler<DataType>::MultipleScalarAssemblyFunction local_asm,
                                                                   const int num_scalars,
                                                                   std::vector< std::vector<DataType> >& vec,
                                                                   typename GlobalAssembler<DataType>::QuadratureSelectionFunction q_select ) const
    {
        // Exactly same implementation as Standard Assembly, so we delegate.
        // TODO(staffan): maybe we should derive from StandardGlobalAssembler instead?
        StandardGlobalAssembler<DataType> assembly;
        assembly.set_quadrature_selection_function ( q_select );
        assembly.assemble_multiple_scalar ( space, local_asm, num_scalars, vec );
    }

    template<class DataType>
    void HpFemAssembler<DataType>::assemble_vector_impl ( const VectorSpace<DataType>& space,
                                                          typename GlobalAssembler<DataType>::VectorAssemblyFunction local_asm,
                                                          typename GlobalAssembler<DataType>::GlobalVector& vec,
                                                          typename GlobalAssembler<DataType>::QuadratureSelectionFunction q_select ) const
    {
        vec.begin_update ( );
        HpVectorAssembly<DataType> assembly ( space, vec );
        assembly.assemble ( local_asm, q_select );
        vec.end_update ( );
    }

    template<class DataType>
    void HpFemAssembler<DataType>::assemble_matrix_impl ( const VectorSpace<DataType>& space,
                                                          typename GlobalAssembler<DataType>::MatrixAssemblyFunction local_asm,
                                                          typename GlobalAssembler<DataType>::GlobalMatrix& mat,
                                                          typename GlobalAssembler<DataType>::QuadratureSelectionFunction q_select ) const
    {
        mat.begin_update ( );
        HpMatrixAssembly<DataType> assembly ( space, mat );
        assembly.assemble ( local_asm, q_select );
        mat.end_update ( );
    }

    template<class DataType>
    void HpFemAssembler<DataType>::assemble_scalar_boundary_impl (
                                                                   const VectorSpace<DataType>& space,
                                                                   typename GlobalAssembler<DataType>::BoundaryScalarAssemblyFunction local_asm,
                                                                   std::vector<DataType>& vec,
                                                                   typename GlobalAssembler<DataType>::FacetQuadratureSelectionFunction fq_select ) const
    {
        // TODO (is this same as StandardAssembler?)
        NOT_YET_IMPLEMENTED;
    }

    template<class DataType>
    void HpFemAssembler<DataType>::assemble_vector_boundary_impl (
                                                                   const VectorSpace<DataType>& space,
                                                                   typename GlobalAssembler<DataType>::BoundaryVectorAssemblyFunction local_asm,
                                                                   typename GlobalAssembler<DataType>::GlobalVector& vec,
                                                                   typename GlobalAssembler<DataType>::FacetQuadratureSelectionFunction fq_select ) const
    {
        // TODO (is this same as StandardAssembler?)
        NOT_YET_IMPLEMENTED;
    }

    template<class DataType>
    void HpFemAssembler<DataType>::assemble_matrix_boundary_impl (
                                                                   const VectorSpace<DataType>& space,
                                                                   typename GlobalAssembler<DataType>::BoundaryMatrixAssemblyFunction local_asm,
                                                                   typename GlobalAssembler<DataType>::GlobalMatrix& mat,
                                                                   typename GlobalAssembler<DataType>::FacetQuadratureSelectionFunction fq_select ) const
    {
        // TODO (is this same as StandardAssembler?)
        NOT_YET_IMPLEMENTED;
    }
    //////////////// End Implementation of HpFemAssembler /////////////////////////////

    template class HpVectorAssembly<double>;
    template class HpVectorAssembly<float>;

    template class HpMatrixAssembly<double>;
    template class HpMatrixAssembly<float>;

    template class HpFemAssembler<double>;
    template class HpFemAssembler<float>;
}

//////////////// OLD CODE ////////////////

#if 0

///
/// \brief Assemble global vector for a VectorSpace.
///
/// \details Assembles a vector
/// \f$b_i = \int_{\Omega}{f(x, \varphi_i)dx}\f$
/// over the domain defined by the mesh
/// associated to a VectorSpace. The integrand is defined through
/// the LocalIntegrator object, whose assemble_local_vector(const
/// Element&, LocalVector&) function should return the locally
/// assembled vector for each element.
///
/// \param[in] space      the VectorSpace for which the assembly is performed
/// \param[in] local_int  functor that performs local vector assembly
/// \param[out] global_vector  the assembled vector \f$b_i\f$
/// \see concept_assembly
///

template<class LocalAssembler, class VectorType>
void assemble_vector ( const VectorSpace& space,
                       LocalAssembler& local_asm,
                       VectorType& global_vector )
{
    std::ofstream octave ( "check_assembly.m", std::ios_base::app );
    octave.precision ( 16 );

    LOG_INFO ( "assembly", "\n=> Start vector assembly" );
    typedef VectorSpace::MeshEntityIterator CellIterator;
    using hiflow::la::MatrixEntry;
    using hiflow::la::MatrixEntryList;

    global_vector.Zeros ( );
    std::vector<int> dof;
    LocalVector lv;

    std::vector<int> traversal_order;
    sort_elements ( space, traversal_order );
    const int num_elements = traversal_order.size ( );

    const FEType* prev_fe_type = 0;
    Quadrature<double> quadrature;

    const DofInterpolation& interp = space.dof ( ).dof_interpolation ( );

    for ( int e = 0; e < num_elements; ++e )
    {
        const int elem_index = traversal_order[e];
        Element elem ( space, elem_index );
        const FEType* fe_type = elem.get_fe_type ( 0 );

        if ( prev_fe_type == 0 || fe_type->get_my_id ( ) != prev_fe_type->get_my_id ( ) )
        {
            // TODO: Now chooses quadrature based on fe type of first variable -- improve this
            choose_quadrature ( *fe_type, quadrature );
            prev_fe_type = fe_type;
        }

        // get global dof indices
        elem.get_dof_indices ( dof );
        const int num_dofs = dof.size ( );

        // assemble locally
        lv.clear ( );
        lv.resize ( num_dofs, 0. );
        local_asm.initialize_for_element ( elem, quadrature );
        local_asm.assemble_local_vector ( elem, lv );

        octave << "% Element " << elem_index << "\n";
        octave << "dof = [" << string_from_range ( dof.begin ( ), dof.end ( ) ) << "] + 1;\n"
                << "b_local = [" << precise_string_from_range ( lv.begin ( ), lv.end ( ) ) << "]';\n"
                << "b(dof) += b_local;\n";

        LOG_INFO ( "assembly", "Element " << elem_index << "\n"
                   << "Dofs = " << string_from_range ( dof.begin ( ), dof.end ( ) ) << "\n"
                   << "lv = " << string_from_range ( lv.begin ( ), lv.end ( ) ) << "\n\n" );

        // assemble into global system
        for ( int i = 0; i < num_dofs; ++i )
        {
            DofInterpolation::const_iterator it = interp.find ( dof[i] );
            if ( it != interp.end ( ) )
            {
                // dof[i] is constrained -> add contributions to dependent dofs
                for ( ConstraintIterator c_it = it->second.begin ( ),
                      c_end = it->second.end ( ); c_it != c_end; ++c_it )
                {
                    if ( space.dof ( ).is_dof_on_sd ( c_it->first ) )
                    {
                        global_vector.Add ( c_it->first, c_it->second * lv[i] );
                    }
                }
            }
            else
            {
                // dof[i] is unconstrained -> add contribution to this dof
                if ( space.dof ( ).is_dof_on_sd ( dof[i] ) )
                {
                    global_vector.Add ( dof[i], lv[i] );
                }
            }
        }
    }
    LOG_INFO ( "assembly", "\n=> End vector assembly" );
    octave << "\n\n";
    octave.close ( );
}

///
/// \brief Assemble global matrix for a VectorSpace.
///
/// \details Assembles a matrix
/// \f$A_{ij} = \int_{\Omega}{f(x,\varphi_i, \varphi_j)dx}\f$
/// over the domain defined by the mesh
/// associated to a VectorSpace. The integrand is defined through
/// the LocalIntegrator object, whose assemble_local_matrix(const
/// Element&, LocalVector&) function should return the locally
/// assembled vector for each element.
///
/// \param[in] space      the VectorSpace for which the assembly is performed
/// \param[in] local_int  functor that performs local vector assembly
/// \param[out] global_matrix  the assembled matrix \f$A_{ij}\f$
/// \see concept_assembly
///

template<class LocalAssembler, class MatrixType>
void assemble_matrix ( const VectorSpace& space,
                       LocalAssembler& local_asm,
                       MatrixType& global_matrix )
{
    std::ofstream octave ( "check_assembly.m", std::ios_base::app );
    octave.precision ( 16 );
    typedef VectorSpace::MeshEntityIterator CellIterator;
    const int dim = space.get_dim ( );

    const DofInterpolation& interp = space.dof ( ).dof_interpolation ( );

    global_matrix.Zeros ( );

    octave << "% ==== Global matrix assembly ====\n";
    octave << "A = zeros(" << global_matrix.nrows_global ( ) << ");\n";

    std::vector<int> dof;
    LocalMatrix lm;

    std::vector<int> traversal_order;
    sort_elements ( space, traversal_order );

    //     std::cout << "Element order = " << string_from_range(traversal_order.begin(), traversal_order.end()) << "\n";

    const int num_elements = traversal_order.size ( );

    // TODO: It would be nice to be able to iterate like this...
    // for (ElementIterator it = space.begin(); it != space.end(); ++it) {
    const FEType* prev_fe_type = 0;
    Quadrature<double> quadrature;

    for ( int e = 0; e < num_elements; ++e )
    {
        const int elem_index = traversal_order[e];
        Element elem ( space, elem_index );

        // TODO: Now chooses quadrature based on fe type of first variable -- improve this
        const FEType* fe_type = elem.get_fe_type ( 0 );

        if ( prev_fe_type == 0 || fe_type->get_my_id ( ) != prev_fe_type->get_my_id ( ) )
        {
            choose_quadrature ( *fe_type, quadrature );
            prev_fe_type = fe_type;
        }

        // get global dof indices
        elem.get_dof_indices ( dof );
        const int num_dofs = dof.size ( );

        // assemble locally
        lm.Resize ( num_dofs, num_dofs );
        lm.Zeros ( );
        local_asm.initialize_for_element ( elem, quadrature );
        local_asm.assemble_local_matrix ( elem, lm );

        octave << "% Element " << elem_index << "\n";
        octave << "dof = [" << string_from_range ( dof.begin ( ), dof.end ( ) ) << "] + 1;\n"
                << "A_local = " << lm << ";\n"
                << "A(dof, dof) += A_local;\n";

        LOG_INFO ( "assembly", "Element " << elem_index << "\n"
                   << "Dofs = " << string_from_range ( dof.begin ( ), dof.end ( ) ) << "\n"
                   << "lm = " << lm << "\n\n" );

        // Assemble into global system.  Only add entries to
        // unconstrained rows and columns, and only if the dof
        // corresponding to the row belongs to the local subdomain.
        for ( int i = 0; i < num_dofs; ++i )
        {
            DofInterpolation::const_iterator it_i = interp.find ( dof[i] );

            if ( it_i != interp.end ( ) )
            {
                // dof[i] is constrained -> add contributions to dependent rows
                for ( ConstraintIterator c_it = it_i->second.begin ( ),
                      c_end = it_i->second.end ( ); c_it != c_end; ++c_it )
                {
                    if ( space.dof ( ).is_dof_on_sd ( c_it->first ) )
                    {

                        for ( int j = 0; j < num_dofs; ++j )
                        {
                            DofInterpolation::const_iterator it_j = interp.find ( dof[j] );

                            if ( it_j != interp.end ( ) )
                            {
                                // dof[j] is constrained -> add attributions to dependent columns
                                // TODO: are these not cleared at the end anyway?
                                for ( ConstraintIterator c2_it = it_j->second.begin ( ),
                                      c2_end = it_j->second.end ( ); c2_it != c2_end; ++c2_it )
                                {
                                    global_matrix.Add ( c_it->first,
                                                        c2_it->first,
                                                        c_it->second * c2_it->second * lm ( i, j ) );
                                }
                            }
                            else
                            {
                                // dof[j] unconstrained -> add contribution to dof[j] column
                                global_matrix.Add ( c_it->first,
                                                    dof[j],
                                                    c_it->second * lm ( i, j ) );
                            }
                        }
                    }
                }
            }
            else
            {
                // dof[i] is unconstrained
                if ( space.dof ( ).is_dof_on_sd ( dof[i] ) )
                {
                    for ( int j = 0; j < num_dofs; ++j )
                    {
                        DofInterpolation::const_iterator it_j = interp.find ( dof[j] );
                        if ( it_j != interp.end ( ) )
                        {
                            for ( ConstraintIterator c_it = it_j->second.begin ( ),
                                  c_end = it_j->second.end ( ); c_it != c_end; ++c_it )
                            {
                                // dof[j] is constrained -> add attributions to dependent columns
                                global_matrix.Add ( dof[i],
                                                    c_it->first,
                                                    c_it->second * lm ( i, j ) );
                            }
                        }
                        else
                        {
                            // dof[j] unconstrained - assemble normally
                            global_matrix.Add ( dof[i], dof[j], lm ( i, j ) );
                        }
                    }
                }
            }
        }
    }

    // Set rows of constrained dofs to identity to obtain non-singular
    // matrix
    std::vector<int> constrained_dofs;
    for ( DofInterpolation::const_iterator it = interp.begin ( ), end = interp.end ( );
          it != end; ++it )
    {

        if ( space.dof ( ).is_dof_on_sd ( it->first ) )
        {
            constrained_dofs.push_back ( it->first );
        }
    }

    {
        const int DEBUG_LEVEL = 3;
        LOG_DEBUG ( 3, "Constrained dofs in assemble_matrix() =\n"
                    << string_from_range ( constrained_dofs.begin ( ), constrained_dofs.end ( ) ); )
    }

    if ( !constrained_dofs.empty ( ) )
    {
        global_matrix.ZeroRows ( &constrained_dofs.front ( ), constrained_dofs.size ( ), 1. );
    }

    octave << "\n\n";
    octave.close ( );
}
#endif
