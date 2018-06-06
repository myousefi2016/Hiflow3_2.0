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
#include "common/pointers.h"
#include "common/sorted_array.h"
#include "dof/dof_interpolation.h"
#include "fem/fetype.h"
#include "linear_algebra/couplings.h"
#include "linear_algebra/coupled_matrix.h"
#include "space/element.h"
#include "space/vector_space.h"
#include "common/log.h"

#include <algorithm>
#include <numeric>

namespace hiflow
{
    using namespace doffem;
    //using doffem::FEType;

    template<class DataType>
    void initialize_linear_algebra_objects ( const VectorSpace<DataType>& space,
                                             const GlobalAssembler<DataType>& global_asm,
                                             la::Couplings<DataType>& couplings,
                                             typename GlobalAssembler<DataType>::GlobalMatrix& matrix )
    {
        SparsityStructure sparsity;
        global_asm.compute_sparsity_structure ( space, sparsity );

        couplings.InitializeCouplings ( sparsity.off_diagonal_rows, sparsity.off_diagonal_cols );

        matrix.InitStructure (
                               vec2ptr ( sparsity.diagonal_rows ),
                               vec2ptr ( sparsity.diagonal_cols ),
                               sparsity.diagonal_rows.size ( ),
                               vec2ptr ( sparsity.off_diagonal_rows ),
                               vec2ptr ( sparsity.off_diagonal_cols ),
                               sparsity.off_diagonal_rows.size ( ) );
    }

    template<class DataType>
    void interpolate_constrained_vector ( const VectorSpace<DataType>& space, typename GlobalAssembler<DataType>::GlobalVector& vector )
    {
        // TODO: necessary to take sub-domain into account here?

        const doffem::DofInterpolation& interp = space.dof ( ).dof_interpolation ( );

        const size_t num_constrained_dofs = interp.size ( );

        // return early if there are no constrained dofs
        if ( num_constrained_dofs == 0 )
        {
            return;
        }

        std::vector<int> constrained_dofs;
        constrained_dofs.reserve ( num_constrained_dofs );
        std::vector<DataType> constrained_values;
        constrained_values.reserve ( num_constrained_dofs );

        std::vector<int> dependencies;
        std::vector<DataType> coefficients;
        std::vector<DataType> dependency_values;

        for ( doffem::DofInterpolation::const_iterator it = interp.begin ( ), end = interp.end ( ); it != end; ++it )
        {
            if ( space.dof ( ).is_dof_on_sd ( it->first ) )
            {

                const size_t num_dependencies = it->second.size ( );
                // probably should not happen, but we check to avoid later problems
                if ( num_dependencies > 0 )
                {

                    dependencies.resize ( num_dependencies, -1 );
                    coefficients.resize ( it->second.size ( ), 0. );
                    dependency_values.resize ( num_dependencies, 0. );

                    for ( size_t i = 0; i != num_dependencies; ++i )
                    {
                        dependencies[i] = it->second[i].first; // id of dependencies
                        coefficients[i] = it->second[i].second; // coefficient of dependencies
                    }

                    // get values of dependency dofs from vector
                    vector.GetValues ( &dependencies.front ( ), num_dependencies, &dependency_values.front ( ) );

                    // compute dot product of dependency_values and coefficients
                    DataType val = 0.;
                    for ( size_t i = 0; i != num_dependencies; ++i )
                    {
                        val += coefficients[i] * dependency_values[i];
                    }

                    // store information
                    constrained_dofs.push_back ( it->first );
                    constrained_values.push_back ( val );
                }
                else
                {
                    LOG_INFO ( "Constrained DoF without dependencies found", true );
                    exit ( -1 );
                }
            }

        }

        // write constrained dofs to vector
        vector.SetValues ( &constrained_dofs.front ( ), constrained_dofs.size ( ), &constrained_values.front ( ) );
    }

    //////////////// Implementation GlobalAssembler ////////////////

    template<class DataType>
    GlobalAssembler<DataType>::GlobalAssembler ( )
    : q_select_ ( default_select_ ), fq_select_ ( default_facet_select_ ), should_reset_assembly_target_ ( true )
    {
    }

    template<class DataType>
    GlobalAssembler<DataType>::~GlobalAssembler ( )
    {
    }

    template<class DataType>
    void GlobalAssembler<DataType>::compute_sparsity_structure ( const VectorSpace<DataType>& space,
                                                                 SparsityStructure& sparsity,
                                                                 std::vector<std::vector<bool> > *coupling_vars ) const
    {
        if ( coupling_vars->empty ( ) )
        {
            coupling_vars->resize ( space.get_nb_var ( ) );
            for ( size_t i = 0, i_e = space.get_nb_var ( ); i != i_e; ++i )
            {
                ( *coupling_vars )[i].resize ( space.get_nb_var ( ), true );
            }
        }

        // Assert correct size of coupling_vars

        assert ( coupling_vars->size ( ) == space.get_nb_var ( ) );
        for ( size_t i = 0, i_e = space.get_nb_var ( ); i != i_e; ++i )
        {
            assert ( ( *coupling_vars )[i].size ( ) == space.get_nb_var ( ) );
        }

        compute_sparsity_structure_impl ( space, sparsity, coupling_vars );
    }

    template<class DataType>
    void GlobalAssembler<DataType>::integrate_scalar ( const VectorSpace<DataType>& space,
                                                       ScalarAssemblyFunction local_asm,
                                                       DataType& integral ) const
    {
        const size_t num_cells = space.mesh ( ).num_entities ( space.mesh ( ).tdim ( ) );
        std::vector<DataType> cell_values ( num_cells, 0. );

        assemble_scalar_impl ( space, local_asm, cell_values, q_select_ );
        integral = std::accumulate ( cell_values.begin ( ), cell_values.end ( ), 0. );
    }

    template<class DataType>
    void GlobalAssembler<DataType>::integrate_multiple_scalar ( const VectorSpace<DataType>& space,
                                                                MultipleScalarAssemblyFunction local_asm, const int num_scalars,
                                                                std::vector<DataType>& integral ) const
    {
        const size_t num_cells = space.mesh ( ).num_entities ( space.mesh ( ).tdim ( ) );
        std::vector< std::vector<DataType> > cell_values ( num_cells );
        for ( int l = 0; l < num_cells; ++l )
            cell_values[l].resize ( num_scalars, 0. );

        assemble_multiple_scalar_impl ( space, local_asm, num_scalars, cell_values, q_select_ );
        integral.resize ( num_scalars, 0. );

        for ( int i = 0; i < num_cells; ++i )
        {
            for ( int l = 0; l < num_scalars; ++l )
            {
                integral[l] += cell_values[i][l];
            }
        }
    }

    template<class DataType>
    void GlobalAssembler<DataType>::assemble_scalar ( const VectorSpace<DataType>& space,
                                                      ScalarAssemblyFunction local_asm,
                                                      std::vector<DataType>& vec ) const
    {
        if ( should_reset_assembly_target_ )
        {
            vec.clear ( );
            const size_t num_cells = space.mesh ( ).num_entities ( space.mesh ( ).tdim ( ) );
            vec.resize ( num_cells, 0. );
        }

        assemble_scalar_impl ( space, local_asm, vec, q_select_ );
    }

    template<class DataType>
    void GlobalAssembler<DataType>::assemble_multiple_scalar ( const VectorSpace<DataType>& space,
                                                               MultipleScalarAssemblyFunction local_asm, const int num_scalars,
                                                               std::vector< std::vector<DataType> >& vec ) const
    {
        if ( should_reset_assembly_target_ )
        {
            vec.clear ( );
            const size_t num_cells = space.mesh ( ).num_entities ( space.mesh ( ).tdim ( ) );
            vec.resize ( num_cells );
            for ( int l = 0; l < num_cells; ++l )
            {
                vec[l].resize ( num_scalars, 0. );
            }
        }

        assemble_multiple_scalar_impl ( space, local_asm, num_scalars, vec, q_select_ );
    }

    template<class DataType>
    void GlobalAssembler<DataType>::assemble_vector ( const VectorSpace<DataType>& space,
                                                      VectorAssemblyFunction local_asm,
                                                      GlobalVector& vec ) const
    {
        if ( should_reset_assembly_target_ )
        {
            vec.Zeros ( );
        }

        assemble_vector_impl ( space, local_asm, vec, q_select_ );
    }

    template<class DataType>
    void GlobalAssembler<DataType>::assemble_matrix ( const VectorSpace<DataType>& space,
                                                      MatrixAssemblyFunction local_asm,
                                                      GlobalMatrix& matrix ) const
    {
        if ( should_reset_assembly_target_ )
        {
            matrix.Zeros ( );
        }

        assemble_matrix_impl ( space, local_asm, matrix, q_select_ );
    }

    template<class DataType>
    void GlobalAssembler<DataType>::integrate_scalar_boundary ( const VectorSpace<DataType>& space,
                                                                BoundaryScalarAssemblyFunction local_asm,
                                                                DataType& integral ) const
    {
        // TODO: actually the number of boundary facets is needed
        const size_t num_facets = space.mesh ( ).num_entities ( space.mesh ( ).tdim ( ) - 1 );
        std::vector<DataType> cell_values ( num_facets, 0. );

        assemble_scalar_boundary_impl ( space, local_asm, cell_values, fq_select_ );
        integral = std::accumulate ( cell_values.begin ( ), cell_values.end ( ), 0. );
    }

    template<class DataType>
    void GlobalAssembler<DataType>::maximize_scalar_boundary ( const VectorSpace<DataType>& space,
                                                               BoundaryScalarAssemblyFunction local_asm,
                                                               DataType& maximum ) const
    {
        // TODO: actually the number of boundary facets is needed
        const size_t num_facets = space.mesh ( ).num_entities ( space.mesh ( ).tdim ( ) - 1 );
        std::vector<DataType> cell_values ( num_facets, 0. );

        assemble_scalar_boundary_impl ( space, local_asm, cell_values, fq_select_ );
        maximum = *std::max_element ( cell_values.begin ( ), cell_values.end ( ) );
    }

    template<class DataType>
    void GlobalAssembler<DataType>::assemble_scalar_boundary ( const VectorSpace<DataType>& space,
                                                               BoundaryScalarAssemblyFunction local_asm,
                                                               std::vector<DataType>& vec ) const
    {
        if ( should_reset_assembly_target_ )
        {
            vec.clear ( );
            const size_t num_cells = space.mesh ( ).num_entities ( space.mesh ( ).tdim ( ) - 1 );
            vec.resize ( num_cells, 0. );
        }

        assemble_scalar_boundary_impl ( space, local_asm, vec, fq_select_ );
    }

    template<class DataType>
    void GlobalAssembler<DataType>::assemble_vector_boundary ( const VectorSpace<DataType>& space,
                                                               BoundaryVectorAssemblyFunction local_asm,
                                                               GlobalVector& vec ) const
    {
        if ( should_reset_assembly_target_ )
        {
            vec.Zeros ( );
        }

        assemble_vector_boundary_impl ( space, local_asm, vec, fq_select_ );
    }

    template<class DataType>
    void GlobalAssembler<DataType>::assemble_matrix_boundary ( const VectorSpace<DataType>& space,
                                                               BoundaryMatrixAssemblyFunction local_asm,
                                                               GlobalMatrix& matrix ) const
    {
        if ( should_reset_assembly_target_ )
        {
            matrix.Zeros ( );
        }

        assemble_matrix_boundary_impl ( space, local_asm, matrix, fq_select_ );
    }

    template<class DataType>
    void GlobalAssembler<DataType>::set_quadrature_selection_function ( QuadratureSelectionFunction q_select )
    {
        if ( q_select != 0 )
        {
            q_select_ = q_select;
        }
        else
        {
            q_select_ = default_select_;
        }
    }

    template<class DataType>
    void GlobalAssembler<DataType>::should_reset_assembly_target ( bool should_reset )
    {
        should_reset_assembly_target_ = should_reset;
    }

    /// The default quadrature selection chooses a quadrature rule that is accurate to 3 * max(fe_degree).

    template<class DataType>
    struct DefaultQuadratureSelection
    {

        DefaultQuadratureSelection ( ) : last_fe_id_ ( FEType<DataType>::NOT_SET ), last_order_ ( -1 )
        {
        }

        void operator() ( const Element<DataType>& elem, Quadrature<DataType>& quadrature )
        {
            const typename FEType<DataType>::FiniteElement fe_id = elem.get_fe_type ( 0 )->get_my_id ( );

            int fe_deg = 0;

            // compute maxmimum FE degree for all variables
            for ( size_t v = 0, end_v = elem.get_num_variables ( ); v != end_v; ++v )
            {
                fe_deg = std::max ( fe_deg, elem.get_fe_type ( v )->get_fe_deg ( ) );
            }

            const int desired_order = 3 * fe_deg;

            // Return early if we already have the desired quadrature.
            // This is a very important optimization, since setting
            // the quadrature is quite expensive. The elements are
            // typically sorted before traversal, which means that we
            // can minimize the number of quadrature switches through this.
            if ( quadrature.size ( ) > 0 && fe_id == last_fe_id_ && desired_order == last_order_ )
            {
                return;
            }

            switch ( fe_id )
            {
                case FEType<DataType>::LAGRANGE_LINE:
                    quadrature.set_cell_type ( 1 );
                    quadrature.set_quadrature_by_order ( "GaussLine", desired_order );
                    break;
                case FEType<DataType>::LAGRANGE_TRI:
                    quadrature.set_cell_type ( 2 );
                    quadrature.set_quadrature_by_order ( "GaussTriangle", desired_order );
                    break;
                case FEType<DataType>::LAGRANGE_QUAD:
                    quadrature.set_cell_type ( 3 );
                    quadrature.set_quadrature_by_order ( "GaussQuadrilateral", desired_order );
                    break;
                case FEType<DataType>::LAGRANGE_TET:
                    quadrature.set_cell_type ( 4 );
                    quadrature.set_quadrature_by_order ( "GaussTetrahedron", desired_order );
                    break;
                case FEType<DataType>::LAGRANGE_HEX:
                    quadrature.set_cell_type ( 5 );
                    quadrature.set_quadrature_by_order ( "GaussHexahedron", desired_order );
                    break;
                case FEType<DataType>::LAGRANGE_PYR:
                    quadrature.set_cell_type ( 6 );
                    quadrature.set_quadrature_by_order ( "GaussPyramid", desired_order );
                    break;
                default:
                    assert ( false );
            };

            last_fe_id_ = fe_id;
            last_order_ = desired_order;
        }

        typename FEType<DataType>::FiniteElement last_fe_id_;
        int last_order_;
    };

    /// The default facet quadrature selection chooses a facet quadrature rule that is accurate to 2 * max(fe_degree).

    template<class DataType>
    struct DefaultFacetQuadratureSelection
    {

        DefaultFacetQuadratureSelection ( ) : last_fe_id_ ( FEType<DataType>::NOT_SET ), last_order_ ( -1 )
        {
        }

        void operator() ( const Element<DataType>& elem, Quadrature<DataType>& quadrature, int facet_number )
        {
            const typename FEType<DataType>::FiniteElement fe_id = elem.get_fe_type ( 0 )->get_my_id ( );

            int fe_deg = 0;

            // compute maxmimum FE degree for all variables
            for ( size_t v = 0, end_v = elem.get_num_variables ( ); v != end_v; ++v )
            {
                fe_deg = std::max ( fe_deg, elem.get_fe_type ( v )->get_fe_deg ( ) );
            }

            const int desired_order = 3 * fe_deg;

            // Only set base quadrature if it changed.
            // This is a very important optimization, since setting
            // the quadrature is quite expensive. The elements are
            // typically sorted before traversal, which means that we
            // can minimize the number of quadrature switches through this.
            if ( !( base_quadrature_.size ( ) > 0 && fe_id == last_fe_id_ && desired_order == last_order_ ) || fe_id == 6 )
            {
                switch ( fe_id )
                {
                    case FEType<DataType>::LAGRANGE_LINE:
                        assert ( 0 );
                        break;
                    case FEType<DataType>::LAGRANGE_TRI:
                        base_quadrature_.set_cell_type ( 2 );
                        base_quadrature_.set_quadrature_by_order ( "GaussLine", desired_order );
                        break;
                    case FEType<DataType>::LAGRANGE_QUAD:
                        base_quadrature_.set_cell_type ( 3 );
                        base_quadrature_.set_quadrature_by_order ( "GaussLine", desired_order );
                        break;
                    case FEType<DataType>::LAGRANGE_TET:
                        base_quadrature_.set_cell_type ( 4 );
                        base_quadrature_.set_quadrature_by_order ( "GaussTriangle", desired_order );
                        break;
                    case FEType<DataType>::LAGRANGE_HEX:
                        base_quadrature_.set_cell_type ( 5 );
                        base_quadrature_.set_quadrature_by_order ( "GaussQuadrilateral", desired_order );
                        break;
                    case FEType<DataType>::LAGRANGE_PYR:
                        base_quadrature_.set_cell_type ( 6 );
                        if ( facet_number == 0 )
                        {
                            base_quadrature_.set_quadrature_by_order ( "GaussQuadrilateral", desired_order );
                        }
                        else
                        {
                            base_quadrature_.set_quadrature_by_order ( "GaussTriangle", desired_order );
                        }
                        break;
                    default:
                        assert ( false );
                };

                last_fe_id_ = fe_id;
                last_order_ = desired_order;
            }
            quadrature.set_facet_quadrature ( base_quadrature_, base_quadrature_.get_cell_type ( ), facet_number );
        }

        typename FEType<DataType>::FiniteElement last_fe_id_;
        int last_order_;
        Quadrature<DataType> base_quadrature_;
    };

    template<>
    const GlobalAssembler<double>::QuadratureSelectionFunction
    GlobalAssembler<double>::default_select_ = DefaultQuadratureSelection<double>( );
    template<>
    const GlobalAssembler<double>::FacetQuadratureSelectionFunction
    GlobalAssembler<double>::default_facet_select_ = DefaultFacetQuadratureSelection<double>( );

    template<>
    const GlobalAssembler<float>::QuadratureSelectionFunction
    GlobalAssembler<float>::default_select_ = DefaultQuadratureSelection<float>( );
    template<>
    const GlobalAssembler<float>::FacetQuadratureSelectionFunction
    GlobalAssembler<float>::default_facet_select_ = DefaultFacetQuadratureSelection<float>( );

    template class GlobalAssembler<double>;
    template class GlobalAssembler<float>;

    template void interpolate_constrained_vector<double>( const VectorSpace<double>& space, GlobalAssembler<double>::GlobalVector& vector );
    template void interpolate_constrained_vector<float>( const VectorSpace<float>& space, GlobalAssembler<float>::GlobalVector& vector );

    /*template void initialize_linear_algebra_objects<double>(const VectorSpace& space, const GlobalAssembler<double>& global_asm,
                                                            la::Couplings& couplings,
                                                            typename GlobalAssembler<double>::GlobalMatrix& matrix);
    template void initialize_linear_algebra_objects<float>(const VectorSpace& space, const GlobalAssembler<float>& global_asm,
                                                            la::Couplings& couplings,
                                                            typename GlobalAssembler<float>::GlobalMatrix& matrix);*/

}
