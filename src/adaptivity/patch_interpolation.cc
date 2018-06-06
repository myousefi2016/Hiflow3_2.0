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

#include "patch_interpolation.h"
#include "common/log.h"
#include "common/vector_algebra.h"
#include "dof/degree_of_freedom.h"
#include "fem/fetype.h"

const int DEBUG_LEVEL = 1;

namespace hiflow
{
    ////////////////////////////////////////////////////
    /////// SpacePatchInterpolation /////////////////////
    ////////////////////////////////////////////////////

    template<class LAD, IMPL MESH_IMPL>
    SpacePatchInterpolation<LAD, MESH_IMPL>::SpacePatchInterpolation ( )
    : eps_ ( 1e-12 ),
    num_var_ ( 0 ),
    initialized_ ( false ),
    mesh_ ( NULL ),
    input_mesh_ ( NULL ),
    input_space_ ( NULL ),
    tdim_ ( 0 ),
    gdim_ ( 0 ),
    rank_ ( 0 )
    {
        this->degrees_.resize ( num_var_ );
        this->is_cg_.resize ( num_var_ );
#ifndef WITH_P4EST
        LOG_ERROR ( "SpacePatchInterpolation: Need P4EST to use!" );
        exit ( -1 );
#endif
    }

    template<class LAD, IMPL MESH_IMPL>
    void SpacePatchInterpolation<LAD, MESH_IMPL>::init ( VectorSpace<DataType> const * input_space )
    {
        this->clear ( );
        this->input_mesh_ = input_space->meshPtr ( );
        this->tdim_ = input_mesh_->tdim ( );
        this->gdim_ = input_mesh_->gdim ( );
        this->input_space_ = input_space;
        this->num_var_ = input_space->get_nb_var ( );

        assert ( this->input_mesh_ != NULL );

        const MPI_Comm& comm = input_space->get_mpi_comm ( );

        MPI_Comm_rank ( comm, &this->rank_ );

        // check if patch interpolation is possible for given space
        bool valid_space = this->check_space ( );
        if ( !valid_space )
        {
            LOG_DEBUG ( 0, " SpacePatchInterpolation: input space is not valid ! " );
            exit ( -1 );
        }

        // get fe degrees of input space
        this->degrees_.resize ( this->num_var_, 0 );
        this->is_cg_ .resize ( this->num_var_, true );

        for ( int v = 0; v<this->num_var_; ++v )
        {
            this->is_cg_[v] = input_space->fe_manager ( ).get_ca ( v );
            this->degrees_[v] = input_space->fe_manager ( ).get_fe_on_cell ( 0, v )->get_fe_deg ( ) * 2;
        }

        // create coarse mesh

        this->copy_mesh ( );

        // coarsen mesh uniformly
        assert ( this->input_mesh_->is_uniformly_coarsenable ( ) );

        std::vector<int> refs ( this->mesh_->num_entities ( this->tdim_ ), -1 );
        this->mesh_ = this->mesh_->refine ( refs );
        this->compute_ghost ( );

        // initialize interpolating space
        this->space_.Init ( this->degrees_, *this->mesh_ );
        this->couplings_.Init ( comm, this->space_.dof ( ) );
        SparsityStructure sparsity;
        this->global_asm_.compute_sparsity_structure ( this->space_, sparsity );
        this->couplings_.InitializeCouplings ( sparsity.off_diagonal_rows, sparsity.off_diagonal_cols );

        // create DOF mapping from interpolating space to input space
        this->create_dof_mapping ( );

        this->initialized_ = true;
    }

    template<class LAD, IMPL MESH_IMPL>
    bool SpacePatchInterpolation<LAD, MESH_IMPL>::check_space ( ) const
    {
        // check if mesh is uniformly coarsenable
        if ( !this->input_mesh_->is_uniformly_coarsenable ( ) )
        {
            LOG_DEBUG ( 0, "Mesh is not unifromly coarsenable " );
            return false;
        }

        // check if fe functions have same degree on every cell
        for ( int v = 0; v<this->num_var_; ++v )
        {
            int ref_deg = this->input_space_->fe_manager ( ).get_fe_on_cell ( 0, v )->get_fe_deg ( );
            for ( EntityNumber jc = 1; jc<this->input_mesh_->num_entities ( this->tdim_ ); ++jc )
            {
                if ( this->input_space_->fe_manager ( ).get_fe_on_cell ( jc, v )->get_fe_deg ( ) != ref_deg )
                {
                    LOG_DEBUG ( 0, "Non-uniform FE degree " );
                    return false;
                }
            }
        }
        return true;
    }

    template<class LAD, IMPL MESH_IMPL>
    void SpacePatchInterpolation<LAD, MESH_IMPL>::create_dof_mapping ( )
    {
        std::map<Id, std::vector<Id> > descendants;

        // Loop over all cells in fine mesh
        for ( EntityIterator cell_it = this->input_mesh_->begin ( this->tdim_ ); cell_it != this->input_mesh_->end ( this->tdim_ ); cell_it++ )
        {
            Id child_id = cell_it->id ( );
            EntityNumber child_index = cell_it->index ( );
            Id parent_id = this->input_mesh_->get_parent_cell_id ( child_index );

            if ( parent_id >= 0 )
            {
                std::map<Id, std::vector<Id> >::iterator it = descendants.find ( parent_id );
                if ( it == descendants.end ( ) )
                {
                    std::vector<Id> tmp_id ( 1, child_id );
                    descendants[parent_id] = tmp_id;
                }
                else
                {
                    it->second.push_back ( child_id );
                }
            }
        }

        // Loop over all cells in coarse mesh <-> higher order space
        for ( EntityIterator cell_it = this->mesh_->begin ( this->tdim_ ); cell_it != this->mesh_->end ( this->tdim_ ); cell_it++ )
        {
            // Get children ids -> corresponding cells in fine mesh
            Id parent_id = cell_it->id ( );
            EntityNumber parent_index = cell_it->index ( );

            if ( !this->mesh_->cell_is_local ( parent_index ) )
            {
                continue;
            }

            std::vector<Id> children_ids = descendants[parent_id];

            // get higher order dofs on coarse cell
            std::vector<int> high_dof_ids;
            this->space_.GetDofIndices ( 0, *cell_it, &high_dof_ids );

            int num_high_dofs = high_dof_ids.size ( );
            std::vector< Coord > high_coords;
            this->space_.dof ( ).get_coord_on_cell ( 0, cell_it->index ( ), high_coords );

            LOG_DEBUG ( 2, " number of children " << children_ids.size ( ) );
            // loop over children cell entities
            for ( int l = 0; l < children_ids.size ( ); ++l )
            {
                int child_index = -1;
                for ( int index = 0; index<this->input_mesh_->num_entities ( this->tdim_ ); ++index )
                {
                    if ( this->input_mesh_->get_id ( this->tdim_, index ) == children_ids[l] )
                    {
                        child_index = index;
                        break;
                    }
                }
                Entity child_cell = this->input_mesh_->get_entity ( this->tdim_, child_index );

                // get lower order dofs on fine cell
                std::vector<int> low_dof_ids;
                this->input_space_->GetDofIndices ( 0, child_cell, &low_dof_ids );
                int num_low_dofs = low_dof_ids.size ( );
                std::vector< Coord > low_coords;
                this->input_space_->dof ( ).get_coord_on_cell ( 0, child_cell.index ( ), low_coords );

                // DOF identification:
                // PERF_TODO: avoid identification by physical coordinates
                // Loop over all dofs in parent_cell
                for ( int i = 0; i < num_high_dofs; i++ )
                {
                    if ( rank_ != this->input_space_->dof ( ).owner_of_dof ( high_dof_ids.at ( i ) ) )
                    {
                        continue;
                    }

                    // Loop over all dofs in current child cell
                    for ( int j = 0; j < num_low_dofs; ++j )
                    {
                        if ( rank_ != this->space_.dof ( ).owner_of_dof ( low_dof_ids.at ( j ) ) )
                        {
                            //continue;
                        }

                        double dist = 0;
                        for ( int d = 0; d<this->gdim_; ++d )
                        {
                            dist += ( high_coords[i][d] - low_coords[j][d] ) * ( high_coords[i][d] - low_coords[j][d] );
                        }
                        dist = std::sqrt ( dist );

                        LOG_DEBUG ( 2, " distance " << dist );
                        // check if high and low dof coincide
                        if ( dist <= this->eps_ )
                        {
                            this->high2low_[high_dof_ids[i]] = low_dof_ids[j];
                            break;
                        }
                    }
                }
            }
        }
    }

    template<class LAD, IMPL MESH_IMPL>
    void SpacePatchInterpolation<LAD, MESH_IMPL>::interpolate ( const VectorType& input_vector, VectorType& vector ) const
    {
        if ( !this->initialized_ )
        {
            std::cout << "SpacePatchInterpolation is not initialized ! " << std::endl;
            exit ( -1 );
        }

        // Loop over all higher order dofs and copy df values from low order vector to higher order vector
        for ( std::map<int, int>::const_iterator it = this->high2low_.begin ( ); it != this->high2low_.end ( ); ++it )
        {
            int high_id = it->first;
            int low_id = it->second;
            double low_val;

            low_val = input_vector.GetValue ( low_id );
            vector.SetValues ( &high_id, 1, &low_val );
        }
        vector.Update ( );
    }

    template<class LAD, IMPL MESH_IMPL>
    void SpacePatchInterpolation<LAD, MESH_IMPL>::copy_mesh ( )
    {
        if ( MESH_IMPL == mesh::IMPL_P4EST )
        {
            this->mesh_ = new MeshPXest ( this->tdim_, this->gdim_ );
        }
        else
        {
            exit ( -1 );
        }
        this->mesh_->deep_copy_from ( this->input_mesh_ );
    }

    template<class LAD, IMPL MESH_IMPL>
    void SpacePatchInterpolation<LAD, MESH_IMPL>::compute_ghost ( )
    {
        SharedVertexTable shared_verts;
        const MPI_Comm& comm = this->input_space_->get_mpi_comm ( );

        this->mesh_ = compute_ghost_cells ( *this->mesh_, comm, shared_verts, MESH_IMPL, this->input_mesh_->get_ghost_layer_width ( ) );
    }

    template<class LAD, IMPL MESH_IMPL>
    void SpacePatchInterpolation<LAD, MESH_IMPL>::clear ( )
    {
        this->num_var_ = 0;
        this->eps_ = 1e-12;
        this->degrees_.clear ( );
        this->is_cg_.clear ( );
        this->high2low_.clear ( );
        this->tdim_ = 0;
        this->gdim_ = 0;
        this->input_mesh_ = NULL;
        this->input_space_ = NULL;
        this->rank_ = 0;

        // clear space
        this->space_.Clear ( );
        this->couplings_.Clear ( );

        // TODO delete mesh?
        this->mesh_ = NULL;

        this->initialized_ = false;
    }

    template class SpacePatchInterpolation<LADescriptorCoupledD, mesh::IMPL_P4EST>;
#ifdef WITH_HYPRE
    template class SpacePatchInterpolation<LADescriptorHypreD, IMPL_P4EST>;
#endif

    ////////////////////////////////////////////////////
    /////// TimePatchInterpolation /////////////////////
    ////////////////////////////////////////////////////

    template<typename T, class DataType>
    T TimePatchInterpolation<T, DataType>::constant ( const DataType c, const T& val_prev, const T& val, const T& val_next ) const
    {
        switch ( this->rel_time_ )
        {
            case 0:
                return val;
                break;
            case 1:
                return val_next;
                break;
        }
    }

    template<typename T, class DataType>
    T TimePatchInterpolation<T, DataType>::linear ( const DataType c, const T& val_prev, const T& val, const T& val_next ) const
    {
        switch ( this->rel_time_ )
        {
            case 0:
                return (c * val + ( 1. - c ) * val_prev );
                break;
            case 1:
                return (c * val_next + ( 1. - c ) * val );
                break;
        }
    }

    template<typename T, class DataType>
    T TimePatchInterpolation<T, DataType>::dt_linear ( const DataType c, const T& val_prev, const T& val, const T& val_next ) const
    {
        switch ( this->rel_time_ )
        {
            case 0:
                return (val - val_prev ) * ( 1. / this->dT_pc_ );
                break;
            case 1:
                return (val_next - val ) * ( 1. / this->dT_cn_ );
                break;
        }
    }

    template<typename T, class DataType>
    T TimePatchInterpolation<T, DataType>::patch_linear ( const DataType c, const T& val_prev, const T& val, const T& val_next ) const
    {
        DataType t = this->get_absolut_time ( c );

        return (this->b_c1_ * t + this->b_c0_ ) * val
                + ( this->b_n1_ * t + this->b_n0_ ) * val_next;
    }

    template<typename T, class DataType>
    T TimePatchInterpolation<T, DataType>::patch_dt_linear ( const DataType c, const T& val_prev, const T& val, const T& val_next ) const
    {
        return (this->b_c1_ * val + this->b_n1_ * val_next );
    }

    template<typename T, class DataType>
    T TimePatchInterpolation<T, DataType>::patch_quadratic ( const DataType c, const T& val_prev, const T& val, const T& val_next ) const
    {
        DataType t = this->get_absolut_time ( c );
        DataType tt = t*t;

        return (this->a_n2_ * tt + this->a_n1_ * t ) * val_next
                + ( this->a_c2_ * tt + this->a_c1_ * t ) * val
                + ( this->a_p2_ * tt + this->a_p1_ * t + this->a_p0_ ) * val_prev;
    }

    template<typename T, class DataType>
    T TimePatchInterpolation<T, DataType>::patch_dt_quadratic ( const DataType c, const T& val_prev, const T& val, const T& val_next ) const
    {
        DataType t = this->get_absolut_time ( c );

        return (2. * this->a_n2_ * t + this->a_n1_ ) * val_next
                + ( 2. * this->a_c2_ * t + this->a_c1_ ) * val
                + ( 2. * this->a_p2_ * t + this->a_p1_ ) * val_prev;
    }

    template<typename T, class DataType>
    void TimePatchInterpolation<T, DataType>::set_time_steps ( DataType dT_pc, DataType dT_cn, int rel_time )
    {
        this->dT_pc_ = dT_pc;
        this->dT_cn_ = dT_cn;
        this->rel_time_ = rel_time;
    }

    template<typename T, class DataType>
    void TimePatchInterpolation<T, DataType>::compute_weight_tau_coeff ( )
    {
        DataType t1 = this->dT_pc_;
        DataType t2 = this->dT_pc_ + this->dT_cn_;
        DataType alpha = t2 * t2 - t1*t2;

        this->a_n2_ = 1. / alpha;
        this->a_n1_ = -t1 / alpha;
        this->a_c2_ = -t2 / ( t1 * alpha );
        this->a_c1_ = t2 * t2 / ( t1 * alpha ); //1. / t1 + t2 / alpha;
        this->a_p2_ = ( t2 - t1 ) / ( t1 * alpha ); //t2 / t1 - 1.;
        this->a_p1_ = -( t1 + t2 ) / ( t1 * t2 ); //-1. / t1 + t1 / alpha - t2 / alpha;
        this->a_p0_ = 1.;

        this->b_n1_ = 2. / t2;
        this->b_n0_ = -t1 / t2;
        this->b_c1_ = -2. / t2;
        this->b_c0_ = 1. + t1 / t2;
    }

    template<typename T, class DataType>
    DataType TimePatchInterpolation<T, DataType>::get_absolut_time ( DataType c ) const
    {
        switch ( this->rel_time_ )
        {
            case 0:
                return this->dT_pc_ * c;
                break;
            case 1:
                return this->dT_pc_ + c * this->dT_cn_;
                break;
        }

    }

    template<typename T, class DataType>
    void TimePatchInterpolation<T, DataType>::clear ( )
    {
        dT_pc_ = 0.;
        dT_cn_ = 0.;
        rel_time_ = 0;

        a_n2_ = 0.;
        a_n1_ = 0.;
        a_c2_ = 0.;
        a_c1_ = 0.;
        a_p2_ = 0.;
        a_p1_ = 0.;
        a_p0_ = 0.;

        b_n1_ = 0.;
        b_n0_ = 0.;
        b_c1_ = 0.;
        b_c0_ = 0.;
    }

    template class TimePatchInterpolation< double, double>;
    template class TimePatchInterpolation< Vec<2, double>, double>;
    template class TimePatchInterpolation< Vec<3, double>, double>;
    template class TimePatchInterpolation< Mat<2, 2, double>, double>;
    template class TimePatchInterpolation< Mat<3, 3, double>, double>;
}
