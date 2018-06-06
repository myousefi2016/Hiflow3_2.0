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

#include "dg_assembly.h"
#include "../linear_algebra/coupled_matrix.h"

#include "mesh/mesh.h"
#include "common/sort_permutation.h"
#include "common/log.h"

const int DEBUG_LEVEL = 0;
const int PRINT_RANK = 0;

/// \author Staffan Ronnas, Jonathan Schwegler, Simon Gawlok

namespace hiflow
{

    using namespace mesh;
    using namespace doffem;

    // Add entries in local matrix to global matrix.

    template<class DataType>
    void add_global ( const VectorSpace<DataType>& space,
                      const std::vector<int>& row_dofs,
                      const std::vector<int>& col_dofs,
                      const typename GlobalAssembler<DataType>::LocalMatrix& lm,
                      typename GlobalAssembler<DataType>::GlobalMatrix& gm )
    {

        const int num_dofs_row = row_dofs.size ( );
        const int num_dofs_col = col_dofs.size ( );

        std::vector<int> row_dofs_sort_permutation;
        std::vector<int> col_dofs_sort_permutation;
        std::vector<int> row_dofs_sorted ( num_dofs_row );
        std::vector<int> col_dofs_sorted ( num_dofs_col );

        // get permutation for sorting dofs
        sortingPermutation ( row_dofs, row_dofs_sort_permutation );
        sortingPermutation ( col_dofs, col_dofs_sort_permutation );

        // fill sorted dof array
        for ( size_t i = 0; i != num_dofs_row; ++i )
        {
            row_dofs_sorted[i] = row_dofs[row_dofs_sort_permutation[i]];
        }
        for ( size_t i = 0; i != num_dofs_col; ++i )
        {
            col_dofs_sorted[i] = col_dofs[col_dofs_sort_permutation[i]];
        }

        // create row array
        std::vector<int> row_indices;
        std::vector<int> row_permutation;
        row_indices.reserve ( num_dofs_row );
        row_permutation.reserve ( num_dofs_row );
        for ( size_t i = 0; i != num_dofs_row; ++i )
        {
            if ( space.dof ( ).is_dof_on_sd ( row_dofs_sorted[i] ) )
            {
                row_indices.push_back ( row_dofs_sorted[i] );
                row_permutation.push_back ( row_dofs_sort_permutation[i] );
            }
        }

        // fill reduced and sorted local matrix
        typename GlobalAssembler<DataType>::LocalMatrix local_mat_sorted_reduced;
        if ( row_indices.size ( ) > 0 && col_dofs_sorted.size ( ) > 0 )
        {
            local_mat_sorted_reduced.Resize ( row_indices.size ( ), col_dofs_sorted.size ( ) );
            for ( size_t i = 0; i != row_indices.size ( ); ++i )
            {
                for ( size_t j = 0; j != col_dofs_sorted.size ( ); ++j )
                {
                    local_mat_sorted_reduced ( i, j ) = lm ( row_permutation[i], col_dofs_sort_permutation[j] );
                }
            }

            // Add local to global matrix
            gm.Add ( vec2ptr ( row_indices ), row_indices.size ( ), vec2ptr ( col_dofs_sorted ), col_dofs_sorted.size ( ), &local_mat_sorted_reduced ( 0, 0 ) );
        }
    }

    template<class DataType>
    void add_global ( const VectorSpace<DataType>& space,
                      const std::vector<int>& dofs,
                      const typename GlobalAssembler<DataType>::LocalVector& lv,
                      typename GlobalAssembler<DataType>::GlobalVector& vec )
    {

        const int num_dofs = dofs.size ( );

        std::vector<int> dofs_sort_permutation;
        std::vector<int> dofs_sorted ( num_dofs );

        // get permutation for sorting dofs
        sortingPermutation ( dofs, dofs_sort_permutation );

        // fill sorted dof array
        for ( size_t i = 0; i != num_dofs; ++i )
        {
            dofs_sorted[i] = dofs[dofs_sort_permutation[i]];
        }

        // create row array
        std::vector<int> row_indices;
        row_indices.reserve ( num_dofs );

        typename GlobalAssembler<DataType>::LocalVector local_vec_sorted;
        local_vec_sorted.reserve ( num_dofs );

        for ( size_t i = 0; i != num_dofs; ++i )
        {
            if ( space.dof ( ).is_dof_on_sd ( dofs_sorted[i] ) )
            {
                row_indices.push_back ( dofs_sorted[i] );
                local_vec_sorted.push_back ( lv[dofs_sort_permutation[i]] );
            }
        }

        // Add local to global vector
        if ( row_indices.size ( ) > 0 )
        {
            vec.Add ( vec2ptr ( row_indices ), row_indices.size ( ), vec2ptr ( local_vec_sorted ) );
        }
    }

    template<class DataType>
    void init_master_quadrature ( const Element<DataType>& slave_elem,
                                  const Element<DataType>& master_elem,
                                  const Quadrature<DataType>& slave_quad,
                                  Quadrature<DataType>& master_quad )
    {
        // Exit early if both elements are the same.
        if ( master_elem.get_cell ( ).id ( ) == slave_elem.get_cell ( ).id ( ) )
        {
            master_quad = slave_quad; // copy quadrature (might not be necessary)
            return;
        }

        const doffem::CellTransformation<DataType>* Ts = slave_elem.get_cell_transformation ( );
        const doffem::CellTransformation<DataType>* Tm = master_elem.get_cell_transformation ( );

        const int num_q = slave_quad.size ( );
        std::vector<DataType> xc ( num_q, 0. ), yc ( num_q, 0. ), zc ( num_q, 0. ), wgt ( num_q, 0. );

        const int dim = slave_elem.get_cell ( ).tdim ( );

        std::vector<DataType> ref_pt ( dim ), phys_pt ( dim );

        if ( dim == 2 )
        {
            for ( size_t q = 0; q != num_q; ++q )
            {
                // reference point
                ref_pt[0] = slave_quad.x ( q );
                ref_pt[1] = slave_quad.y ( q );

                // physical point
                phys_pt[0] = Ts->x ( ref_pt );
                phys_pt[1] = Ts->y ( ref_pt );

                // reference point on master cell
                Tm->inverse ( phys_pt[0], phys_pt[1], xc[q], yc[q] );

                // weight
                wgt[q] = slave_quad.w ( q );
            }

        }
        else if ( dim == 3 )
        {
            for ( size_t q = 0; q != num_q; ++q )
            {
                // reference point

                ref_pt[0] = slave_quad.x ( q );
                ref_pt[1] = slave_quad.y ( q );
                ref_pt[2] = slave_quad.z ( q );

                // physical point
                phys_pt[0] = Ts->x ( ref_pt );
                phys_pt[1] = Ts->y ( ref_pt );
                phys_pt[2] = Ts->z ( ref_pt );
                // reference point on master cell
                Tm->inverse ( phys_pt[0], phys_pt[1], phys_pt[2], xc[q], yc[q], zc[q] );

                // weight
                wgt[q] = slave_quad.w ( q );
            }
        }
        master_quad.set_custom_quadrature ( xc, yc, zc, wgt );
    }

    template <class DataType>
    struct DefaultInterfaceQuadratureSelection
    {

        void operator() ( const Element<DataType>& master_elem,
                const Element<DataType>& slave_elem,
                int master_facet_number,
                int slave_facet_number,
                Quadrature<DataType>& master_quad,
                Quadrature<DataType>& slave_quad )
        {

            const typename FEType<DataType>::FiniteElement fe_id = slave_elem.get_fe_type ( 0 )->get_my_id ( );

            int fe_deg = 0;

            // compute maxmimum FE degree for all variables
            for ( size_t v = 0, end_v = slave_elem.get_num_variables ( ); v != end_v; ++v )
            {
                int var_fe_deg = slave_elem.get_fe_type ( v )->get_fe_deg ( );
                fe_deg = std::max ( fe_deg, var_fe_deg );
            }
            for ( size_t v = 0, end_v = master_elem.get_num_variables ( ); v != end_v; ++v )
            {
                int var_fe_deg = master_elem.get_fe_type ( v )->get_fe_deg ( );
                fe_deg = std::max ( fe_deg, var_fe_deg );
            }

            const int desired_order = 3 * fe_deg;
            if ( !( base_quad_.size ( ) > 0 && fe_id == last_fe_id_ && desired_order == last_order_ ) )
            {
                switch ( fe_id )
                {
                    case FEType<DataType>::LAGRANGE_TRI:
                        base_quad_.set_cell_type ( 1 );
                        base_quad_.set_quadrature_by_order ( "GaussLine", desired_order );
                        break;
                    case FEType<DataType>::LAGRANGE_QUAD:
                        base_quad_.set_cell_type ( 1 );
                        base_quad_.set_quadrature_by_order ( "GaussLine", desired_order );
                        break;
                    case FEType<DataType>::LAGRANGE_TET:
                        base_quad_.set_cell_type ( 2 );
                        base_quad_.set_quadrature_by_order ( "GaussTriangle", desired_order );
                        break;
                    case FEType<DataType>::LAGRANGE_HEX:
                        base_quad_.set_cell_type ( 3 );
                        base_quad_.set_quadrature_by_order ( "EconomicalGaussQuadrilateral", desired_order );
                        break;
                    case FEType<DataType>::LAGRANGE_PYR:
                        if ( slave_facet_number == 0 )
                        {
                            base_quad_.set_cell_type ( 3 );
                            base_quad_.set_quadrature_by_order ( "EconomicalGaussQuadrilateral", desired_order );
                        }
                        else
                        {
                            base_quad_.set_cell_type ( 2 );
                            base_quad_.set_quadrature_by_order ( "GaussTriangle", desired_order );
                        }
                        break;
                    default:
                        assert ( false );
                };
                last_fe_id_ = fe_id;
                last_order_ = desired_order;
            }
            slave_quad.set_facet_quadrature ( base_quad_,
                                              slave_elem.get_cell ( ).cell_type ( ).tag ( ),
                                              slave_facet_number );
            init_master_quadrature ( slave_elem, master_elem, slave_quad, master_quad );

        }
        Quadrature<DataType> base_quad_;
        typename FEType<DataType>::FiniteElement last_fe_id_;
        int last_order_;

    };

    template <class DataType>
    DGGlobalAssembler<DataType>::DGGlobalAssembler ( )
    : if_q_select_ ( DefaultInterfaceQuadratureSelection<DataType>( ) )
    {
        // By default, we don't reset the target.
        this->should_reset_assembly_target ( false );
    }

    template <class DataType>
    void DGGlobalAssembler<DataType>::assemble_interface_matrix ( const VectorSpace<DataType>& space,
                                                                  InterfaceMatrixAssemblyFun local_asm,
                                                                  typename GlobalAssembler<DataType>::GlobalMatrix& matrix ) const
    {
        if ( this->should_reset_assembly_target_ )
        {
            matrix.Zeros ( );
        }
        // Create interface list from mesh
        ConstMeshPtr mesh = &space.mesh ( ); // OK since pointer is intrusive
        InterfaceList if_list = InterfaceList::create ( mesh );

        typename GlobalAssembler<DataType>::LocalMatrix L_MM, L_MS, L_SM, L_SS;

        std::vector<int> master_dofs, slave_dofs;

        Quadrature<DataType> master_quadrature, slave_quadrature;

        // Loop over interfaces

        for ( InterfaceList::const_iterator it = if_list.begin ( ),
              end_it = if_list.end ( ); it != end_it; ++it )
        {
            int remote_index_master = -10;
            mesh->get_attribute_value ( "_remote_index_", mesh->tdim ( ), it->master_index ( ), &remote_index_master );

            // Master dofs
            Element<DataType> master_elem ( space, it->master_index ( ) );

            space.GetDofIndices ( master_elem.get_cell ( ), &master_dofs );
            const int master_facet_number = it->master_facet_number ( );

            L_MM.Resize ( master_dofs.size ( ), master_dofs.size ( ) );
            L_MM.Zeros ( );

            // Initialize master quadrature
            if_q_select_ ( master_elem, master_elem,
                           master_facet_number, master_facet_number,
                           master_quadrature, master_quadrature );

            const int num_slaves = it->num_slaves ( );
            if ( remote_index_master == -1 )
            {
                if ( num_slaves > 0 )
                {
                    local_asm ( master_elem, master_elem,
                                master_quadrature, master_quadrature,
                                master_facet_number, master_facet_number,
                                INTERFACE_MASTER, INTERFACE_MASTER, L_MM );
                }
                else
                {
                    // boundary facet
                    local_asm ( master_elem, master_elem,
                                master_quadrature, master_quadrature,
                                master_facet_number, master_facet_number,
                                INTERFACE_MASTER, INTERFACE_BOUNDARY, L_MM );
                }
                add_global ( space, master_dofs, master_dofs, L_MM, matrix );
            }
            // Loop over slaves
            for ( size_t s = 0; s != num_slaves; ++s )
            {
                int remote_index_slave = -10;
                mesh->get_attribute_value ( "_remote_index_", mesh->tdim ( ), it->slave_index ( s ), &remote_index_slave );
                Element<DataType> slave_elem ( space, it->slave_index ( s ) );
                space.GetDofIndices ( slave_elem.get_cell ( ), &slave_dofs );
                const int slave_facet_number = it->slave_facet_number ( s );

                // Initialize slave quadrature. NB: only called once per slave.
                if_q_select_ ( master_elem, slave_elem,
                               master_facet_number, slave_facet_number,
                               master_quadrature, slave_quadrature );

                if ( remote_index_master == -1 )
                {
                    // master / slave
                    L_MS.Resize ( master_dofs.size ( ), slave_dofs.size ( ) );
                    L_MS.Zeros ( );
                    local_asm ( slave_elem, master_elem,
                                slave_quadrature, master_quadrature,
                                slave_facet_number, master_facet_number,
                                INTERFACE_SLAVE, INTERFACE_MASTER, L_MS );
                    add_global ( space, master_dofs, slave_dofs, L_MS, matrix );
                }

                if ( remote_index_slave == -1 )
                {
                    // slave / master
                    L_SM.Resize ( slave_dofs.size ( ), master_dofs.size ( ) );
                    L_SM.Zeros ( );
                    local_asm ( master_elem, slave_elem,
                                master_quadrature, slave_quadrature,
                                master_facet_number, slave_facet_number,
                                INTERFACE_MASTER, INTERFACE_SLAVE, L_SM );
                    add_global ( space, slave_dofs, master_dofs, L_SM, matrix );

                    // slave / slave
                    L_SS.Resize ( slave_dofs.size ( ), slave_dofs.size ( ) );
                    L_SS.Zeros ( );
                    local_asm ( slave_elem, slave_elem,
                                slave_quadrature, slave_quadrature,
                                slave_facet_number, slave_facet_number,
                                INTERFACE_SLAVE, INTERFACE_SLAVE, L_SS );
                    add_global ( space, slave_dofs, slave_dofs, L_SS, matrix );
                }
            }
        }
    }

    template <class DataType>
    void DGGlobalAssembler<DataType>::assemble_interface_vector ( const VectorSpace<DataType>& space,
                                                                  InterfaceVectorAssemblyFun local_asm,
                                                                  typename GlobalAssembler<DataType>::GlobalVector& vec ) const
    {
        if ( this->should_reset_assembly_target_ )
        {
            vec.Zeros ( );
        }
        // Create interface list from mesh
        ConstMeshPtr mesh = &space.mesh ( ); // OK since pointer is intrusive

        InterfaceList if_list = InterfaceList::create ( mesh );

        typename GlobalAssembler<DataType>::LocalVector L_MM, L_MS, L_SM, L_SS;

        std::vector<int> master_dofs, slave_dofs;

        Quadrature<DataType> master_quadrature, slave_quadrature;

        // Loop over interfaces

        for ( InterfaceList::const_iterator it = if_list.begin ( ),
              end_it = if_list.end ( ); it != end_it; ++it )
        {
            int remote_index_master = -10;
            mesh->get_attribute_value ( "_remote_index_", mesh->tdim ( ), it->master_index ( ), &remote_index_master );

            // Master dofs
            Element<DataType> master_elem ( space, it->master_index ( ) );

            space.GetDofIndices ( master_elem.get_cell ( ), &master_dofs );
            const int master_facet_number = it->master_facet_number ( );

            L_MM.clear ( );
            L_MM.resize ( master_dofs.size ( ) );

            // Initialize master quadrature
            if_q_select_ ( master_elem, master_elem,
                           master_facet_number, master_facet_number,
                           master_quadrature, master_quadrature );

            const int num_slaves = it->num_slaves ( );
            if ( remote_index_master == -1 )
            {
                if ( num_slaves > 0 )
                {
                    local_asm ( master_elem, master_elem,
                                master_quadrature, master_quadrature,
                                master_facet_number, master_facet_number,
                                INTERFACE_MASTER, INTERFACE_MASTER, L_MM );
                }
                else
                {
                    // boundary facet
                    local_asm ( master_elem, master_elem,
                                master_quadrature, master_quadrature,
                                master_facet_number, master_facet_number,
                                INTERFACE_MASTER, INTERFACE_BOUNDARY, L_MM );
                }
                add_global ( space, master_dofs, L_MM, vec );
            }
            // Loop over slaves
            for ( size_t s = 0; s != num_slaves; ++s )
            {
                int remote_index_slave = -10;
                mesh->get_attribute_value ( "_remote_index_", mesh->tdim ( ), it->slave_index ( s ), &remote_index_slave );
                Element<DataType> slave_elem ( space, it->slave_index ( s ) );
                space.GetDofIndices ( slave_elem.get_cell ( ), &slave_dofs );
                const int slave_facet_number = it->slave_facet_number ( s );

                // Initialize slave quadrature. NB: only called once per slave.
                if_q_select_ ( master_elem, slave_elem,
                               master_facet_number, slave_facet_number,
                               master_quadrature, slave_quadrature );

                if ( remote_index_master == -1 )
                {
                    // master / slave
                    L_MS.clear ( );
                    L_MS.resize ( master_dofs.size ( ) );
                    local_asm ( slave_elem, master_elem,
                                slave_quadrature, master_quadrature,
                                slave_facet_number, master_facet_number,
                                INTERFACE_SLAVE, INTERFACE_MASTER, L_MS );
                    add_global ( space, master_dofs, L_MS, vec );
                }

                if ( remote_index_slave == -1 )
                {
                    // slave / master
                    L_SM.clear ( );
                    L_SM.resize ( slave_dofs.size ( ) );
                    local_asm ( master_elem, slave_elem,
                                master_quadrature, slave_quadrature,
                                master_facet_number, slave_facet_number,
                                INTERFACE_MASTER, INTERFACE_SLAVE, L_SM );
                    add_global ( space, slave_dofs, L_SM, vec );

                    // slave / slave
                    L_SS.clear ( );
                    L_SS.resize ( slave_dofs.size ( ) );
                    local_asm ( slave_elem, slave_elem,
                                slave_quadrature, slave_quadrature,
                                slave_facet_number, slave_facet_number,
                                INTERFACE_SLAVE, INTERFACE_SLAVE, L_SS );
                    add_global ( space, slave_dofs, L_SS, vec );
                }
            }
        }
    }

    template<class DataType>
    void DGGlobalAssembler<DataType>::assemble_interface_scalar ( const VectorSpace<DataType>& space,
                                                                  InterfaceScalarAssemblyFun local_asm,
                                                                  std::vector<DataType>& values ) const
    {

        // Create interface list from mesh
        ConstMeshPtr mesh = &space.mesh ( ); // OK since pointer is intrusive

        InterfaceList if_list = InterfaceList::create ( mesh );

        // Clear and create values data structure
        values.clear ( );
        values.resize ( if_list.size ( ), 0. );

        DataType L_MM, L_MS, L_SM, L_SS;

        Quadrature<DataType> master_quadrature, slave_quadrature;

        // Loop over interfaces
        int i = 0;
        for ( InterfaceList::const_iterator it = if_list.begin ( ),
              end_it = if_list.end ( ); it != end_it; ++it )
        {
            int remote_index_master = -10;
            mesh->get_attribute_value ( "_remote_index_", mesh->tdim ( ), it->master_index ( ), &remote_index_master );

            // Master dofs
            Element<DataType> master_elem ( space, it->master_index ( ) );

            const int master_facet_number = it->master_facet_number ( );

            L_MM = 0.;

            // Initialize master quadrature
            if_q_select_ ( master_elem, master_elem,
                           master_facet_number, master_facet_number,
                           master_quadrature, master_quadrature );

            const int num_slaves = it->num_slaves ( );
            if ( remote_index_master == -1 )
            {
                if ( num_slaves > 0 )
                {
                    local_asm ( master_elem, master_elem,
                                master_quadrature, master_quadrature,
                                master_facet_number, master_facet_number,
                                INTERFACE_MASTER, INTERFACE_MASTER, L_MM );
                }
                else
                {
                    // boundary facet
                    local_asm ( master_elem, master_elem,
                                master_quadrature, master_quadrature,
                                master_facet_number, master_facet_number,
                                INTERFACE_MASTER, INTERFACE_BOUNDARY, L_MM );
                }
                values[i] += L_MM;
            }
            // Loop over slaves
            for ( size_t s = 0; s != num_slaves; ++s )
            {
                int remote_index_slave = -10;
                mesh->get_attribute_value ( "_remote_index_", mesh->tdim ( ), it->slave_index ( s ), &remote_index_slave );
                Element<DataType> slave_elem ( space, it->slave_index ( s ) );
                const int slave_facet_number = it->slave_facet_number ( s );

                // Initialize slave quadrature. NB: only called once per slave.
                if_q_select_ ( master_elem, slave_elem,
                               master_facet_number, slave_facet_number,
                               master_quadrature, slave_quadrature );

                if ( remote_index_master == -1 )
                {
                    // master / slave
                    L_MS = 0.;
                    local_asm ( slave_elem, master_elem,
                                slave_quadrature, master_quadrature,
                                slave_facet_number, master_facet_number,
                                INTERFACE_SLAVE, INTERFACE_MASTER, L_MS );
                    values[i] += L_MS;
                }

                if ( remote_index_slave == -1 )
                {
                    // slave / master
                    L_SM = 0.;
                    local_asm ( master_elem, slave_elem,
                                master_quadrature, slave_quadrature,
                                master_facet_number, slave_facet_number,
                                INTERFACE_MASTER, INTERFACE_SLAVE, L_SM );
                    values[i] += L_SM;

                    // slave / slave
                    L_SS = 0.;
                    local_asm ( slave_elem, slave_elem,
                                slave_quadrature, slave_quadrature,
                                slave_facet_number, slave_facet_number,
                                INTERFACE_SLAVE, INTERFACE_SLAVE, L_SS );
                    values[i] += L_SS;
                }
            }
            ++i;
        }
    }

    template<class DataType>
    void DGGlobalAssembler<DataType>::assemble_interface_scalar_cells ( const VectorSpace<DataType>& space,
                                                                        InterfaceScalarAssemblyFun local_asm,
                                                                        std::vector<DataType>& values ) const
    {

        // Create interface list from mesh
        ConstMeshPtr mesh = &space.mesh ( ); // OK since pointer is intrusive

        InterfaceList if_list = InterfaceList::create ( mesh );

        // Clear and create values data structure
        values.clear ( );
        values.resize ( mesh->num_entities ( mesh->tdim ( ) ), 0. );

        DataType L_MM, L_MS, L_SM, L_SS;

        int rank = space.dof ( ).get_my_subdomain ( );
        int print_rank = 0;

        // Loop over interfaces
        for ( InterfaceList::const_iterator it = if_list.begin ( ),
              end_it = if_list.end ( ); it != end_it; ++it )
        {
            int remote_index_master = -10;
            mesh->get_attribute_value ( "_remote_index_", mesh->tdim ( ), it->master_index ( ), &remote_index_master );

            // Master dofs
            Element<DataType> master_elem ( space, it->master_index ( ) );

            const int master_facet_number = it->master_facet_number ( );

            DataType L_MM, L_MS, L_SM, L_SS;

            Quadrature<DataType> master_master_quadrature;

            L_MM = 0.;

            // Initialize master quadrature
            this->if_q_select_ ( master_elem, master_elem,
                                 master_facet_number, master_facet_number,
                                 master_master_quadrature, master_master_quadrature );

            const size_t num_slaves = it->num_slaves ( );

            if ( remote_index_master == -1 )
            {
                if ( num_slaves > 0 )
                {
                    local_asm ( master_elem, master_elem,
                                master_master_quadrature, master_master_quadrature,
                                master_facet_number, master_facet_number,
                                INTERFACE_MASTER, INTERFACE_MASTER, L_MM );
                }
                else
                {
                    // boundary facet
                    local_asm ( master_elem, master_elem,
                                master_master_quadrature, master_master_quadrature,
                                master_facet_number, master_facet_number,
                                INTERFACE_MASTER, INTERFACE_BOUNDARY, L_MM );
                }

                LOG_DEBUG ( 3, "[" << rank << "] Master index: " << it->master_index ( ) << " with remote index " << remote_index_master
                            << ", add to master_cell L_MM=" << L_MM );

                values[it->master_index ( )] += L_MM;
            }

            // Loop over slaves
            for ( size_t s = 0; s != num_slaves; ++s )
            {
                int remote_index_slave = -10;
                mesh->get_attribute_value ( "_remote_index_", mesh->tdim ( ), it->slave_index ( s ), &remote_index_slave );
                Element<DataType> slave_elem ( space, it->slave_index ( s ) );
                const int slave_facet_number = it->slave_facet_number ( s );

                Quadrature<DataType> master_quadrature, slave_quadrature;

                // Initialize slave quadrature. NB: only called once per slave.
                this->if_q_select_ ( master_elem, slave_elem,
                                     master_facet_number, slave_facet_number,
                                     master_quadrature, slave_quadrature );

                if ( remote_index_master == -1 )
                {
                    // master / slave
                    L_MS = 0.;
                    local_asm ( slave_elem, master_elem,
                                slave_quadrature, master_quadrature,
                                slave_facet_number, master_facet_number,
                                INTERFACE_SLAVE, INTERFACE_MASTER, L_MS );

                    if ( num_slaves > 1 && rank == PRINT_RANK )
                    {
                        LOG_DEBUG ( 2, "[" << rank << "] Master index: " << it->master_index ( ) << " with remote index " << remote_index_master
                                    << " and slave index " << it->slave_index ( s ) << " with remote index " << remote_index_slave
                                    << ", add to master_cell L_MS=" << L_MS );
                    }
                    values[it->master_index ( )] += L_MS;
                }

                if ( remote_index_slave == -1 )
                {
                    // slave / master
                    L_SM = 0.;

                    local_asm ( master_elem, slave_elem,
                                master_quadrature, slave_quadrature,
                                master_facet_number, slave_facet_number,
                                INTERFACE_MASTER, INTERFACE_SLAVE, L_SM );

                    if ( num_slaves > 1 && rank == PRINT_RANK )
                    {
                        LOG_DEBUG ( 2, "[" << rank << "] Master index: " << it->master_index ( ) << " with remote index " << remote_index_master
                                    << " and slave index " << it->slave_index ( s ) << " with remote index " << remote_index_slave
                                    << ", add to slave_cell L_SM=" << L_SM );
                    }
                    values[it->slave_index ( s )] += L_SM;

                    // slave / slave
                    L_SS = 0.;
                    local_asm ( slave_elem, slave_elem,
                                slave_quadrature, slave_quadrature,
                                slave_facet_number, slave_facet_number,
                                INTERFACE_SLAVE, INTERFACE_SLAVE, L_SS );

                    LOG_DEBUG ( 3, "[" << rank << "] Master index: " << it->master_index ( ) << " with remote index " << remote_index_master
                                << " and slave index " << it->slave_index ( s ) << " with remote index " << remote_index_slave
                                << ", add to slvae_cell L_SS=" << L_SS );

                    values[it->slave_index ( s )] += L_SS;
                }
            }

        }
    }

    template<class DataType>
    void DGGlobalAssembler<DataType>::assemble_interface_multiple_scalar_cells ( const VectorSpace<DataType>& space,
                                                                                 InterfaceMultipleScalarAssemblyFun local_asm,
                                                                                 const int num_scalars,
                                                                                 std::vector<typename GlobalAssembler<DataType>::LocalVector>& values ) const
    {

        // Create interface list from mesh
        ConstMeshPtr mesh = &space.mesh ( ); // OK since pointer is intrusive

        InterfaceList if_list = InterfaceList::create ( mesh );

        // Clear and create values data structure
        values.clear ( );
        values.resize ( mesh->num_entities ( mesh->tdim ( ) ) );
        for ( int j = 0; j < values.size ( ); ++j )
        {
            values[j].resize ( num_scalars, 0. );
        }

        int rank = space.dof ( ).get_my_subdomain ( );
        int print_rank = 0;

        // Loop over interfaces
        for ( InterfaceList::const_iterator it = if_list.begin ( ),
              end_it = if_list.end ( ); it != end_it; ++it )
        {
            int remote_index_master = -10;
            mesh->get_attribute_value ( "_remote_index_", mesh->tdim ( ), it->master_index ( ), &remote_index_master );

            // Master dofs
            Element<DataType> master_elem ( space, it->master_index ( ) );

            const int master_facet_number = it->master_facet_number ( );

            typename GlobalAssembler<DataType>::LocalVector L_MM, L_MS;

            Quadrature<DataType> master_master_quadrature;

            // Initialize master quadrature
            this->if_q_select_ ( master_elem, master_elem,
                                 master_facet_number, master_facet_number,
                                 master_master_quadrature, master_master_quadrature );

            const size_t num_slaves = it->num_slaves ( );

            // Boundary integral
            if ( remote_index_master == -1 )
            {
                if ( num_slaves == 0 )
                {
                    // boundary facet
                    local_asm ( master_elem, master_elem,
                                master_master_quadrature, master_master_quadrature,
                                master_facet_number, master_facet_number,
                                INTERFACE_MASTER, INTERFACE_BOUNDARY, L_MM );

                    LOG_DEBUG ( 3, "[" << rank << "] Master index: " << it->master_index ( ) << " with remote index " << remote_index_master
                                << ", add to master_cell L_MM=" << string_from_range ( L_MM.begin ( ), L_MM.end ( ) ) );

                    assert ( values[it->master_index ( )].size ( ) == L_MM.size ( ) );
                    for ( int l = 0; l < L_MM.size ( ); ++l )
                    {
                        values[it->master_index ( )][l] += L_MM[l];
                    }
                }
            }

            // Interface integrals
            // Loop over slaves
            for ( size_t s = 0; s != num_slaves; ++s )
            {
                int remote_index_slave = -10;
                mesh->get_attribute_value ( "_remote_index_", mesh->tdim ( ), it->slave_index ( s ), &remote_index_slave );
                Element<DataType> slave_elem ( space, it->slave_index ( s ) );
                const int slave_facet_number = it->slave_facet_number ( s );

                Quadrature<DataType> master_quadrature, slave_quadrature;

                // Initialize slave quadrature. NB: only called once per slave.
                this->if_q_select_ ( master_elem, slave_elem,
                                     master_facet_number, slave_facet_number,
                                     master_quadrature, slave_quadrature );

                if ( remote_index_master == -1 || remote_index_slave == -1 )
                {
                    local_asm ( slave_elem, master_elem,
                                slave_quadrature, master_quadrature,
                                slave_facet_number, master_facet_number,
                                INTERFACE_SLAVE, INTERFACE_MASTER, L_MS );

                    assert ( values[it->master_index ( )].size ( ) == L_MS.size ( ) );
                    assert ( values[it->slave_index ( s )].size ( ) == L_MS.size ( ) );

                    if ( remote_index_master == -1 )
                    {
                        for ( int l = 0; l < L_MS.size ( ); ++l )
                        {
                            values[it->master_index ( )][l] += L_MS[l];
                        }
                    }
                    if ( remote_index_slave == -1 )
                    {
                        for ( int l = 0; l < L_MS.size ( ); ++l )
                        {
                            values[it->slave_index ( s )][l] += L_MS[l];
                        }
                    }
                }
            }
        }
    }

    template<class DataType>
    void DGGlobalAssembler<DataType>::distribute_interface_to_cell_values_naive ( const VectorSpace<DataType>& space,
                                                                                  std::vector<DataType> &cell_values,
                                                                                  const std::vector<DataType> &interface_values ) const
    {
        // Create interface list from mesh
        ConstMeshPtr mesh = &space.mesh ( ); // OK since pointer is intrusive

        // Check compatibility of mesh and cell_values vector
        assert ( mesh->num_entities ( mesh->tdim ( ) ) == cell_values.size ( ) );

        InterfaceList if_list = InterfaceList::create ( mesh );

        // Check compatibility of interface list and interface_values vector
        assert ( if_list.size ( ) == interface_values.size ( ) );

        // Loop over interfaces
        int i = 0;
        for ( InterfaceList::const_iterator it = if_list.begin ( ),
              end_it = if_list.end ( ); it != end_it; ++it )
        {
            int remote_index_master = -10;
            mesh->get_attribute_value ( "_remote_index_", mesh->tdim ( ), it->master_index ( ), &remote_index_master );

            const int num_slaves = it->num_slaves ( );
            if ( remote_index_master == -1 )
            {
                if ( num_slaves > 0 )
                {
                    // Master only gets half of the contribution of interface value
                    cell_values[it->master_index ( )] += 0.5 * interface_values[i];
                }
                else
                {
                    // boundary facet
                    // Master only gets contribution of interface value
                    cell_values[it->master_index ( )] += interface_values[i];
                }
            }

            // weight per slave
            DataType weight_slave = 0.5;
            if ( num_slaves > 0 )
            {
                weight_slave /= num_slaves;
            }
            // Loop over slaves
            for ( size_t s = 0; s != num_slaves; ++s )
            {
                int remote_index_slave = -10;
                mesh->get_attribute_value ( "_remote_index_", mesh->tdim ( ), it->slave_index ( s ), &remote_index_slave );

                if ( remote_index_slave == -1 )
                {
                    cell_values[it->slave_index ( s )] += weight_slave * interface_values[i];
                }
            }
            ++i;
        }
    }

    template <class DataType>
    void DGGlobalAssembler<DataType>::set_interface_quadrature_selection_fun (
                                                                               IFQuadratureSelectionFun q_select )
    {
        if_q_select_ = q_select;
    }

    template <class DataType>
    void DGGlobalAssembler<DataType>::compute_sparsity_structure_impl ( const VectorSpace<DataType>& space,
                                                                        SparsityStructure& sparsity,
                                                                        std::vector<std::vector<bool> > *coupling_vars ) const
    {

        // Create interface list from mesh
        ConstMeshPtr mesh = &space.mesh ( ); // OK since pointer is intrusive
        InterfaceList if_list = InterfaceList::create ( mesh );

        const int ndof_total = space.dof ( ).ndofs_on_sd ( space.dof ( ).get_my_subdomain ( ) );

        // a set of columns for every row
        std::vector< SortedArray<int> > raw_diag ( ndof_total );
        std::vector< SortedArray<int> > raw_offdiag ( ndof_total );

        std::vector<int> dof_list, slave_dofs;

        for ( InterfaceList::const_iterator it = if_list.begin ( ),
              end_it = if_list.end ( ); it != end_it; ++it )
        {
            for ( size_t test_var = 0, e_test_var = space.get_nb_var ( ); test_var != e_test_var; ++test_var )
            {
                // get dof indices
                Element<DataType> master_elem ( space, it->master_index ( ) );
                space.GetDofIndices ( test_var, master_elem.get_cell ( ), &dof_list );

                // loop over trial variables
                for ( size_t trial_var = 0, e_trial_var = space.get_nb_var ( ); trial_var != e_trial_var; ++trial_var )
                {
                    // check whether test_var and trial_var couple
                    if ( ( *coupling_vars )[test_var][trial_var] )
                    {
                        //from master cell
                        if ( trial_var != test_var )
                        {
                            space.GetDofIndices ( trial_var, master_elem.get_cell ( ), &slave_dofs );
                            dof_list.insert ( dof_list.end ( ), slave_dofs.begin ( ), slave_dofs.end ( ) );
                        }
                        //from slave cell
                        for ( size_t s = 0, s_e = it->num_slaves ( ); s != s_e; ++s )
                        {
                            Element<DataType> slave_elem ( space, it->slave_index ( s ) );
                            space.GetDofIndices ( trial_var, slave_elem.get_cell ( ), &slave_dofs );
                            dof_list.insert ( dof_list.end ( ), slave_dofs.begin ( ), slave_dofs.end ( ) );
                        }
                    }
                }
                // All these dofs now potentially couple with one another.
                const int nd = dof_list.size ( );

                for ( size_t i = 0; i != nd; ++i )
                { // rows
                    if ( !space.dof ( ).is_dof_on_sd ( dof_list[i] ) ) continue; // skip remote rows.

                    // get local row dof index
                    int local_dof_i;
                    space.dof ( ).global2local ( dof_list[i], &local_dof_i );

                    for ( size_t j = 0; j != nd; ++j )
                    { // cols
                        // diagonal coupling (my col)
                        if ( space.dof ( ).is_dof_on_sd ( dof_list[j] ) )
                        {
                            // add if coupling is new
                            raw_diag[local_dof_i].find_insert ( dof_list[j] );
                        }
                        else
                        {
                            raw_offdiag[local_dof_i].find_insert ( dof_list[j] );
                        }
                    }
                } // for (int i=0;...
            }
        }

        // post-process to SparsityStructure
        // compute nnz for diagonal and offdiagonal blocks
        int nnz_diag = 0, nnz_offdiag = 0;
        for ( size_t k = 0; k != ndof_total; ++k )
        {
            nnz_diag += raw_diag[k].size ( );
            nnz_offdiag += raw_offdiag[k].size ( );
        }

        sparsity.diagonal_rows.resize ( nnz_diag );
        sparsity.diagonal_cols.resize ( nnz_diag );
        sparsity.off_diagonal_rows.resize ( nnz_offdiag );
        sparsity.off_diagonal_cols.resize ( nnz_offdiag );

        int global_dof_i;

        // copy diagonal sparsity structure
        int k = 0;
        for ( size_t r = 0; r != ndof_total; ++r )
        {
            space.dof ( ).local2global ( r, &global_dof_i );

            for ( SortedArray<int>::const_iterator it = raw_diag[r].begin ( ),
                  end = raw_diag[r].end ( ); it != end; ++it )
            {
                sparsity.diagonal_rows[k] = global_dof_i;
                sparsity.diagonal_cols[k] = *it;
                ++k;
            }
        }
        assert ( k == nnz_diag );

        // copy off-diagonal sparsity structure
        k = 0;
        for ( size_t r = 0; r != ndof_total; ++r )
        {
            space.dof ( ).local2global ( r, &global_dof_i );

            for ( SortedArray<int>::const_iterator it = raw_offdiag[r].begin ( ),
                  end = raw_offdiag[r].end ( ); it != end; ++it )
            {
                sparsity.off_diagonal_rows[k] = global_dof_i;
                sparsity.off_diagonal_cols[k] = *it;
                ++k;
            }
        }
        assert ( k == nnz_offdiag );
    }

    template class DGGlobalAssembler<float>;
    template class DGGlobalAssembler<double>;
}
