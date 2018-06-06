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

#include "dynamic_mesh.h"

#include "mesh/mesh_pXest.h"
#include "mesh/mesh_db_view.h"
#include "mesh/mesh_tools.h"
#include "mesh/writer.h"
#include "common/log.h"
#include "common/sort_permutation.h"
#include <iterator>

const int DEBUG_LEVEL = 1;

namespace hiflow
{

    template<class LAD, IMPL MESH_IMPL, int DIM>
    DynamicMeshHandler<LAD, MESH_IMPL, DIM>::DynamicMeshHandler ( DynamicMeshProblem<LAD, MESH_IMPL, DIM>* problem, MPI_Comm comm )
    : adapt_type_ ( "None" ),
    problem_ ( problem ),
    mesh_ ( NULL ),
    initial_mesh_ ( NULL ),
    time_mesh_ ( NULL ),
    space_ ( NULL ),
    space_dual_ ( NULL ),
    tdim_ ( DIM ),
    gdim_ ( DIM ),
    active_mesh_index_ ( -1 ),
    comm_ ( comm )
    {
        MPI_Comm_rank ( this->comm_, &rank_ );
    }

    template<class LAD, IMPL MESH_IMPL, int DIM>
    void DynamicMeshHandler<LAD, MESH_IMPL, DIM>::clear ( )
    {
        this->vec_primal_primal_.clear ( );
        this->vec_dual_primal_.clear ( );
        this->vec_dual_dual_.clear ( );
        this->vec_est_primal_.clear ( );
        this->vec_est_dual_.clear ( );

        this->mesh_change_times_.clear ( );
        this->mesh_change_steps_.clear ( );

        this->mesh_ = NULL;

        this->time_mesh_ = NULL;
        this->initial_mesh_ = NULL;

        this->tdim_ = DIM;
        this->gdim_ = DIM;
        this->active_mesh_index_ = -1;
        this->adapt_type_ = "None";
        this->space_ = NULL;
        this->space_dual_ = NULL;
        this->mesh_list_.clear ( );

        for ( int l = 0; l< this->space_list_.size ( ); ++l )
        {
            if ( this->space_list_[l] != NULL )
            {
                delete this->space_list_[l];
            }
        }
        this->space_list_.clear ( );

        for ( int l = 0; l< this->space_list_dual_.size ( ); ++l )
        {
            if ( this->space_list_dual_[l] != NULL )
            {
                delete this->space_list_dual_[l];
            }
        }
        this->space_list_dual_.clear ( );

        /*
        for (int l=0; l< this->space_list_tmp_.size(); ++l)
        {
            if (this->space_list_tmp_[l] != NULL)
            {
                delete this->space_list_tmp_[l];
            }
        }*/
        this->space_list_tmp_.clear ( );

        for ( int l = 0; l< this->couplings_list_.size ( ); ++l )
        {
            if ( this->couplings_list_[l] != NULL )
            {
                delete this->couplings_list_[l];
            }
        }
        this->couplings_list_.clear ( );

        for ( int l = 0; l< this->couplings_list_dual_.size ( ); ++l )
        {
            if ( this->couplings_list_dual_[l] != NULL )
            {
                delete this->couplings_list_dual_[l];
            }
        }
        this->couplings_list_dual_.clear ( );

        for ( int l = 0; l<this->fe_interpolator_.size ( ); ++l )
        {
            for ( int k = 0; k<this->fe_interpolator_[l].size ( ); ++k )
            {
                if ( this->fe_interpolator_[l][k] != NULL )
                {
                    delete this->fe_interpolator_[l][k];
                }
            }
        }
        this->fe_interpolator_.clear ( );

        this->coupling_vars_list_.clear ( );
        this->coupling_vars_list_dual_.clear ( );
    }

    template<class LAD, IMPL MESH_IMPL, int DIM>
    void DynamicMeshHandler<LAD, MESH_IMPL, DIM>::set_update_vectors ( int mode, int type, std::vector<VectorType*> vectors )
    {
        if ( mode == 1 )
        {
            if ( type == 1 )
            {
                this->vec_primal_primal_ = vectors;
            }
        }
        else if ( mode == -1 )
        {
            if ( type == 1 )
            {
                this->vec_dual_primal_ = vectors;
            }
            else if ( type == -1 )
            {
                this->vec_dual_dual_ = vectors;
            }
        }
        else if ( mode == 0 )
        {
            if ( type == 1 )
            {
                this->vec_est_primal_ = vectors;
            }
            else if ( type == -1 )
            {
                this->vec_est_dual_ = vectors;
            }
        }
    }

    template<class LAD, IMPL MESH_IMPL, int DIM>
    int DynamicMeshHandler<LAD, MESH_IMPL, DIM>::mesh_index ( DataType time ) const
    {
        if ( time < 0. )
        {
            return -1;
        }
        if ( time > this->time_mesh_->end ( ) + 1e-10 )
        {
            return -1;
        }

        int num_points = this->mesh_change_times_.size ( );
        int cur_point = 0;
        for ( int t = 0; t < num_points; ++t )
        {
            if ( t == num_points - 1 )
            {
                if ( this->mesh_change_times_[t] <= time )
                {
                    cur_point = t + 1;
                    break;
                }
            }
            else
            {
                if ( time < this->mesh_change_times_.data ( )[t + 1] && time >= this->mesh_change_times_.data ( )[t] )
                {
                    cur_point = t + 1;
                    break;
                }
            }
        }
        return cur_point;
    }

    template<class LAD, IMPL MESH_IMPL, int DIM>
    int DynamicMeshHandler<LAD, MESH_IMPL, DIM>::mesh_index ( int time_step ) const
    {
        if ( time_step < 0 )
        {
            return -1;
        }
        if ( time_step > this->time_mesh_->num_intervals ( ) )
        {
            return -1;
        }
        return this->mesh_index ( this->time_mesh_->time ( time_step ) );
    }

    template<class LAD, IMPL MESH_IMPL, int DIM>
    int DynamicMeshHandler<LAD, MESH_IMPL, DIM>::first_step_for_mesh ( int mesh_index ) const
    {
        assert ( mesh_index >= 0 );
        assert ( mesh_index < this->num_mesh ( ) );

        if ( mesh_index == 0 )
        {
            return 0;
        }

        DataType mesh_start_time = this->mesh_change_times_[mesh_index - 1];

        for ( int t = 0; t<this->time_mesh_->num_intervals ( ); ++t )
        {
            if ( this->time_mesh_->time ( t ) >= mesh_start_time )
            {
                return t;
            }
        }
        return -1;
    }

    template<class LAD, IMPL MESH_IMPL, int DIM>
    int DynamicMeshHandler<LAD, MESH_IMPL, DIM>::last_step_for_mesh ( int mesh_index ) const
    {
        assert ( mesh_index >= 0 );
        assert ( mesh_index < this->num_mesh ( ) );

        if ( mesh_index == this->num_mesh ( ) - 1 )
        {
            return this->time_mesh_->num_intervals ( ) - 1;
        }

        return this->first_step_for_mesh ( mesh_index + 1 ) - 1;
    }

    template<class LAD, IMPL MESH_IMPL, int DIM>
    int DynamicMeshHandler<LAD, MESH_IMPL, DIM>::num_cells_by_step ( int time_step ) const
    {
        int mesh_index = this->mesh_index ( time_step );
        return this->num_cells_by_index ( mesh_index );
    }

    template<class LAD, IMPL MESH_IMPL, int DIM>
    int DynamicMeshHandler<LAD, MESH_IMPL, DIM>::num_cells_by_index ( int mesh_index ) const
    {
        if ( mesh_index < 0 )
        {
            return 0;
        }
        if ( mesh_index >= this->num_mesh ( ) )
        {
            return 0;
        }

        return this->mesh_list_[mesh_index]->num_entities ( DIM );
    }

    template<class LAD, IMPL MESH_IMPL, int DIM>
    VectorSpace<typename LAD::DataType>* DynamicMeshHandler<LAD, MESH_IMPL, DIM>::get_space_by_step ( int time_step, int mode ) const
    {
        if ( time_step < 0 )
        {
            return NULL;
        }
        if ( time_step > this->time_mesh_->num_intervals ( ) )
        {
            return NULL;
        }

        int mesh_index = this->mesh_index ( time_step );

        assert ( mesh_index >= 0 );
        assert ( mesh_index < this->mesh_list_.size ( ) );

        if ( mode == 1 )
        {
            return this->space_list_[mesh_index];
        }
        else if ( mode == -1 )
        {
            return this->space_list_dual_[mesh_index];
        }
        else if ( mode == 2 )
        {
            return this->space_list_tmp_[mesh_index];
        }
        return NULL;
    }

    template<class LAD, IMPL MESH_IMPL, int DIM>
    VectorSpace<typename LAD::DataType>* DynamicMeshHandler<LAD, MESH_IMPL, DIM>::get_space_by_index ( int mesh_index, int mode ) const
    {
        assert ( mesh_index >= 0 );
        assert ( mesh_index < this->num_mesh ( ) );

        if ( mode == 1 )
        {
            return this->space_list_[mesh_index];
        }
        else if ( mode == -1 )
        {
            return this->space_list_dual_[mesh_index];
        }
        else if ( mode == 2 )
        {
            return this->space_list_tmp_[mesh_index];
        }
        return NULL;
    }

    template<class LAD, IMPL MESH_IMPL, int DIM>
    Couplings<typename LAD::DataType>* DynamicMeshHandler<LAD, MESH_IMPL, DIM>::get_couplings_by_index ( int mesh_index, int mode ) const
    {
        assert ( mesh_index >= 0 );
        assert ( mesh_index < this->num_mesh ( ) );

        if ( mode == 1 )
        {
            return this->couplings_list_[mesh_index];
        }
        else if ( mode == -1 )
        {
            return this->couplings_list_dual_[mesh_index];
        }
        return NULL;
    }

    template<class LAD, IMPL MESH_IMPL, int DIM>
    MeshPtr DynamicMeshHandler<LAD, MESH_IMPL, DIM>::get_mesh_by_step ( int time_step ) const
    {
        if ( time_step < 0 )
        {
            return NULL;
        }
        if ( time_step > this->time_mesh_->num_intervals ( ) )
        {
            return NULL;
        }
        int mesh_index = this->mesh_index ( time_step );

        return this->get_mesh_by_index ( mesh_index );
    }

    template<class LAD, IMPL MESH_IMPL, int DIM>
    MeshPtr DynamicMeshHandler<LAD, MESH_IMPL, DIM>::get_mesh_by_index ( int mesh_index ) const
    {
        assert ( mesh_index >= 0 );
        assert ( mesh_index < this->num_mesh ( ) );

        return this->mesh_list_[mesh_index];
    }

    template<class LAD, IMPL MESH_IMPL, int DIM>
    bool DynamicMeshHandler<LAD, MESH_IMPL, DIM>::need_to_change_mesh ( DataType last_time, DataType cur_time ) const
    {
        int cur_mesh_index = this->mesh_index ( cur_time );
        int last_mesh_index = this->mesh_index ( last_time );

        if ( cur_mesh_index == last_mesh_index )
        {
            return false;
        }
        return true;
    }

    template<class LAD, IMPL MESH_IMPL, int DIM>
    bool DynamicMeshHandler<LAD, MESH_IMPL, DIM>::need_to_change_mesh ( int last_step, int cur_step ) const
    {
        if ( last_step < 0 )
        {
            return true;
        }

        DataType cur_time = this->time_mesh_->time ( cur_step );
        DataType last_time = this->time_mesh_->time ( last_step );

        return this->need_to_change_mesh ( last_time, cur_time );
    }

    template<class LAD, IMPL MESH_IMPL, int DIM>
    bool DynamicMeshHandler<LAD, MESH_IMPL, DIM>::need_to_change_mesh ( int step ) const
    {
        int mesh_index = this->mesh_index ( step );
        if ( mesh_index != this->active_mesh_index_ )
        {
            return true;
        }
        return false;
    }

    template<class LAD, IMPL MESH_IMPL, int DIM>
    void DynamicMeshHandler<LAD, MESH_IMPL, DIM>::update_mesh_change_steps ( )
    {
        this->mesh_change_steps_.clear ( );
        int old_mesh = -1;
        int new_mesh = -1;
        for ( int t = 0; t<this->time_mesh_->num_intervals ( ); ++t )
        {
            old_mesh = new_mesh;
            new_mesh = this->mesh_index ( t );
            if ( new_mesh >= 0 )
            {
                if ( old_mesh != new_mesh )
                {
                    this->mesh_change_steps_.insert ( t );
                }
            }
        }
    }

    template<class LAD, IMPL MESH_IMPL, int DIM>
    void DynamicMeshHandler<LAD, MESH_IMPL, DIM>::fill_mesh_list ( )
    {
        int num_mesh = this->mesh_change_times_.size ( ) + 1;
        this->mesh_list_.clear ( );
        this->mesh_list_.resize ( num_mesh );

        for ( int l = 0; l < num_mesh; ++l )
        {
            if ( MESH_IMPL == mesh::IMPL_P4EST )
            {
                this->mesh_list_[l] = new MeshPXest ( this->tdim_, this->gdim_ );
            }
            else if ( MESH_IMPL == mesh::IMPL_DBVIEW )
            {
                this->mesh_list_[l] = new MeshDbView ( this->tdim_, this->gdim_ );
            }
            this->mesh_list_[l] ->deep_copy_from ( this->initial_mesh_ );
        }
    }

    template<class LAD, IMPL MESH_IMPL, int DIM>
    void DynamicMeshHandler<LAD, MESH_IMPL, DIM>::init_fe_spaces ( )
    {
        int num_mesh = this->mesh_list_.size ( );
        LOG_DEBUG ( 1, "Init FE spaces for " << num_mesh << " different meshes" );

        int num_space = this->space_list_.size ( );
        int num_space_dual = this->space_list_dual_.size ( );

        this->space_list_tmp_.clear ( );
        for ( int l = 0; l < num_space; ++l )
        {
            if ( this->space_list_[l] != NULL )
            {
                delete this->space_list_[l];
            }
        }

        for ( int l = 0; l < num_space_dual; ++l )
        {
            if ( this->space_list_dual_[l] != NULL )
            {
                delete this->space_list_dual_[l];
            }
        }

        for ( int l = 0; l<this->couplings_list_.size ( ); ++l )
        {
            if ( this->couplings_list_[l] != NULL )
            {
                delete this->couplings_list_[l];
            }
        }

        for ( int l = 0; l<this->couplings_list_dual_.size ( ); ++l )
        {
            if ( this->couplings_list_dual_[l] != NULL )
            {
                delete this->couplings_list_dual_[l];
            }
        }

        this->space_list_ .clear ( );
        this->space_list_dual_ .clear ( );
        this->coupling_vars_list_ .clear ( );
        this->coupling_vars_list_dual_.clear ( );
        this->couplings_list_ .clear ( );
        this->couplings_list_dual_ .clear ( );

        this->space_list_ .resize ( num_mesh );
        this->space_list_dual_ .resize ( num_mesh );
        this->coupling_vars_list_ .resize ( num_mesh );
        this->coupling_vars_list_dual_.resize ( num_mesh );
        this->couplings_list_ .resize ( num_mesh );
        this->couplings_list_dual_ .resize ( num_mesh );

        for ( int l = 0; l < num_mesh; ++l )
        {
            this->space_list_[l] = new VectorSpace<DataType>( );
            this->space_list_dual_[l] = new VectorSpace<DataType>( );

            this->couplings_list_[l] = new Couplings<DataType> ( );
            this->couplings_list_dual_[l] = new Couplings<DataType> ( );

            this->problem_->setup_space ( *this->space_list_[l], this->coupling_vars_list_[l], this->mesh_list_[l], 1 );
            this->problem_->setup_space ( *this->space_list_dual_[l], this->coupling_vars_list_dual_[l], this->mesh_list_[l], -1 );
        }

        this->space_list_tmp_ = this->space_list_;
        this->create_fe_interpolator ( num_mesh );

        LOG_DEBUG ( 1, "Mesh change times " << string_from_range ( this->mesh_change_times_.data ( ).begin ( ), this->mesh_change_times_.data ( ).end ( ) ) );
        for ( int m = 0; m<this->time_mesh_->num_levels ( ); ++m )
        {
            //LOG_DEBUG(1, "Time mesh " << m << " has points: " << string_from_range(this->time_mesh_->get_all_times(m).begin(), this->time_mesh_->get_all_times(m).end()));
        }
    }

    template<class LAD, IMPL MESH_IMPL, int DIM>
    bool DynamicMeshHandler<LAD, MESH_IMPL, DIM>::update ( int time_step, int mode, bool initial_call )
    {
        if ( !this->need_to_change_mesh ( time_step ) && !initial_call )
        {
            LOG_DEBUG ( 1, "No Mesh change needed for time step " << time_step );
            return false;
        }

        // mode: 1 = primal, -1 = dual, 0 = estimation
        int pp_mesh_index = this->mesh_index ( time_step - 2 );
        int p_mesh_index = this->mesh_index ( time_step - 1 );
        int n_mesh_index = this->mesh_index ( time_step + 1 );
        int c_mesh_index = this->mesh_index ( time_step );

        int old_mesh_index = this->active_mesh_index_;
        int new_mesh_index = c_mesh_index;

        LOG_DEBUG ( 1, "Change from mesh " << this->active_mesh_index_ << " to mesh " << new_mesh_index );

        // 1. save previous and (in case of dual problem ) next solution
        std::vector<VectorType*> p_vectors;
        std::vector<VectorType*> d_vectors;

        if ( !initial_call )
        {
            if ( mode == 1 )
            {
                p_vectors.resize ( this->vec_primal_primal_.size ( ) );
                for ( int l = 0; l < p_vectors.size ( ); ++l )
                {
                    p_vectors[l] = new VectorType;
                    p_vectors[l]->CloneFrom ( *this->vec_primal_primal_[l] );
                }
            }
            if ( mode == -1 )
            {
                p_vectors.resize ( this->vec_dual_primal_.size ( ) );
                for ( int l = 0; l < p_vectors.size ( ); ++l )
                {
                    p_vectors[l] = new VectorType;
                    p_vectors[l]->CloneFrom ( *this->vec_dual_primal_[l] );
                }

                d_vectors.resize ( this->vec_dual_dual_.size ( ) );
                for ( int l = 0; l < d_vectors.size ( ); ++l )
                {
                    d_vectors[l] = new VectorType;
                    d_vectors[l]->CloneFrom ( *this->vec_dual_dual_[l] );
                }
            }
            if ( mode == 0 )
            {
                p_vectors.resize ( this->vec_est_primal_.size ( ) );
                for ( int l = 0; l < p_vectors.size ( ); ++l )
                {
                    p_vectors[l] = new VectorType;
                    p_vectors[l]->CloneFrom ( *this->vec_est_primal_[l] );
                }

                d_vectors.resize ( this->vec_est_dual_.size ( ) );
                for ( int l = 0; l < d_vectors.size ( ); ++l )
                {
                    d_vectors[l] = new VectorType;
                    d_vectors[l]->CloneFrom ( *this->vec_est_dual_[l] );
                }
            }
        }

        // 2. change pointers
        this->mesh_ = this->get_mesh_by_index ( c_mesh_index );
        this->space_ = this->get_space_by_index ( c_mesh_index, 1 );
        this->space_dual_ = this->get_space_by_index ( c_mesh_index, -1 );
        this->active_mesh_index_ = c_mesh_index;
        int num_global_cells = this->space_->meshPtr ( )->num_global_cells ( this->space_->get_mpi_comm ( ) );
        this->problem_ ->set_active_space ( this->space_, 1 );
        this->problem_ ->set_active_space ( this->space_dual_, -1 );
        this->problem_ ->set_active_mesh ( this->mesh_ );
        this->problem_ ->set_active_mesh_index ( c_mesh_index );

        LOG_DEBUG ( 2, "  Active space has " << this->space_->dof ( ).ndofs_global ( ) << " total dofs and " << num_global_cells << " cells " );

        // 3. reinit LA objects
        if ( mode == 1 )
        {
            this->problem_->setup_LA_primal ( *this->space_, *this->space_dual_,
                                              this->coupling_vars_list_[c_mesh_index], this->coupling_vars_list_dual_[c_mesh_index],
                                              *this->couplings_list_[c_mesh_index], *this->couplings_list_dual_[c_mesh_index] );
        }
        else if ( mode == -1 )
        {
            this->problem_->setup_LA_dual ( *this->space_, *this->space_dual_,
                                            this->coupling_vars_list_[c_mesh_index], this->coupling_vars_list_dual_[c_mesh_index],
                                            *this->couplings_list_[c_mesh_index], *this->couplings_list_dual_[c_mesh_index] );
        }
        else if ( mode == 0 )
        {
            this->problem_->setup_LA_est ( *this->space_, *this->space_dual_,
                                           this->coupling_vars_list_[c_mesh_index], this->coupling_vars_list_dual_[c_mesh_index],
                                           *this->couplings_list_[c_mesh_index], *this->couplings_list_dual_[c_mesh_index] );
        }

        if ( initial_call )
        {
            return true;
        }

        // 4. interpolate old vectors w.r.t. new space
        if ( mode == 1 )
        {
            for ( int l = 0; l < p_vectors.size ( ); ++l )
            {
                LOG_DEBUG ( 2, "  Interpolate primal vector " << l << " from mesh " << old_mesh_index << " to mesh " << new_mesh_index );
                this->interpolate_vector ( *p_vectors[l], old_mesh_index, *this->vec_primal_primal_[l], new_mesh_index );
            }
        }
        else if ( mode == -1 )
        {
            for ( int l = 0; l < p_vectors.size ( ); ++l )
            {
                LOG_DEBUG ( 2, "  Interpolate primal vector " << l << " from mesh " << old_mesh_index << " to mesh " << new_mesh_index );
                this->interpolate_vector ( *p_vectors[l], old_mesh_index, *this->vec_dual_primal_[l], new_mesh_index );
            }
            for ( int l = 0; l < d_vectors.size ( ); ++l )
            {
                LOG_DEBUG ( 2, "  Interpolate dual vector " << l << " from mesh " << old_mesh_index << " to mesh " << new_mesh_index );
                this->interpolate_vector ( *d_vectors[l], old_mesh_index, *this->vec_dual_dual_[l], new_mesh_index );
            }
        }
        else if ( mode == 0 )
        {
            for ( int l = 0; l < p_vectors.size ( ); ++l )
            {
                LOG_DEBUG ( 2, "  Interpolate primal vector " << l << " from mesh " << old_mesh_index << " to mesh " << new_mesh_index );
                this->interpolate_vector ( *p_vectors[l], old_mesh_index, *this->vec_est_primal_[l], new_mesh_index );
            }
            for ( int l = 0; l < d_vectors.size ( ); ++l )
            {
                LOG_DEBUG ( 2, "  Interpolate dual vector " << l << " from mesh " << old_mesh_index << " to mesh " << new_mesh_index );
                this->interpolate_vector ( *d_vectors[l], old_mesh_index, *this->vec_est_dual_[l], new_mesh_index );
            }
            /*this->prepare_higher_order_space();*/
        }

        // 5. reinit linear solver
        /*
        if (mode == 1)
        {
            if (time_step == 0)
            {
                this->prepare_linear_solver(true);
            }
            else
            {
                this->prepare_linear_solver(false);
            }
        }
        if (mode == -1)
        {
            this->prepare_linear_solver_dual();
        }
         */

        for ( int l = 0; l < p_vectors.size ( ); ++l )
        {
            delete p_vectors[l];
        }
        for ( int l = 0; l < d_vectors.size ( ); ++l )
        {
            delete d_vectors[l];
        }
        return true;
    }

    template<class LAD, IMPL MESH_IMPL, int DIM>
    void DynamicMeshHandler<LAD, MESH_IMPL, DIM>::create_fe_interpolator ( int num_mesh )
    {
        this->fe_interpolator_.resize ( num_mesh );
        for ( int l = 0; l < num_mesh; ++l )
        {
            this->fe_interpolator_[l].resize ( num_mesh );
            for ( int k = 0; k < num_mesh; ++k )
            {
                if ( this->fe_interpolator_[l][k] != NULL )
                {
                    this->fe_interpolator_[l][k]->clear ( );
                }
                else
                {
                    this->fe_interpolator_[l][k] = new FEInterpolationLagrange<LAD, DIM>( );
                }
            }
        }
    }

    template<class LAD, IMPL MESH_IMPL, int DIM>
    void DynamicMeshHandler<LAD, MESH_IMPL, DIM>::interpolate_vector ( const VectorType& in_vec, int in_index, VectorType& out_vec, int out_index )
    {
        if ( !this->fe_interpolator_.at ( in_index ).at ( out_index )->is_initialized ( ) )
        {
            this->fe_interpolator_.at ( in_index ).at ( out_index )->init ( this->get_space_by_index ( in_index, 1 ), this->get_space_by_index ( out_index, 1 ) );
        }

        this->fe_interpolator_.at ( in_index ).at ( out_index )->interpolate ( in_vec, out_vec );
    }

    template<class LAD, IMPL MESH_IMPL, int DIM>
    void DynamicMeshHandler<LAD, MESH_IMPL, DIM>::add_mesh_change_times ( const std::vector<DataType>& additional_changes, std::vector<int>& indicator_mesh_indices )
    {
        if ( this->adapt_type_ == "None" || this->adapt_type_ == "Fixed" )
        {
            this->space_list_tmp_ = this->space_list_;
            return;
        }

        int num_steps = this->time_mesh_->num_intervals ( );
        int num_mesh = this->num_mesh ( );

        SortedArray<DataType> new_change_list;
        new_change_list.data ( ) = this->mesh_change_times_.data ( );
        std::map<DataType, bool> is_new_change;
        for ( int t = 0; t < new_change_list.size ( ); ++t )
        {
            is_new_change.insert ( std::pair<DataType, bool> ( new_change_list.data ( ).at ( t ), false ) );
        }

        // build new mesh change list
        for ( int t = 0; t < additional_changes.size ( ); ++t )
        {
            new_change_list.insert ( additional_changes[t] );
            is_new_change.insert ( std::pair<DataType, bool> ( additional_changes[t], true ) );
        }

        LOG_DEBUG ( 2, "  Old mesh change list: " << string_from_range ( this->mesh_change_times_.begin ( ), this->mesh_change_times_.end ( ) ) );

        // fill in copy of meshes
        std::vector<MeshPtr> new_mesh_list;

        std::vector< SortedArray<int> > old2new_mesh_indices;
        old2new_mesh_indices.resize ( this->num_mesh ( ) );
        old2new_mesh_indices[0].insert ( 0 );

        new_mesh_list.push_back ( this->get_mesh_by_index ( 0 ) );
        space_list_tmp_.clear ( );
        space_list_tmp_.push_back ( this->get_space_by_index ( 0, 1 ) );

        for ( int t = 0; t < new_change_list.size ( ); ++t )
        {
            DataType time = new_change_list[t];
            bool new_mesh = is_new_change[time];
            int mesh_index = this->mesh_index ( time );
            old2new_mesh_indices[mesh_index].insert ( t + 1 );

            LOG_DEBUG ( 2, " mesh change time " << time << " is new mesh " << new_mesh << " gets mesh with index " << mesh_index );

            if ( new_mesh )
            {
                MeshPtr tmp_mesh;
                if ( MESH_IMPL == mesh::IMPL_P4EST )
                {
                    tmp_mesh = new MeshPXest ( DIM, DIM );
                }
                else if ( MESH_IMPL == mesh::IMPL_DBVIEW )
                {
                    tmp_mesh = new MeshDbView ( DIM, DIM );
                }

                tmp_mesh->deep_copy_from ( this->get_mesh_by_index ( mesh_index ) );
                new_mesh_list.push_back ( tmp_mesh );
            }
            else
            {
                new_mesh_list.push_back ( this->get_mesh_by_index ( mesh_index ) );
            }
            space_list_tmp_.push_back ( this->get_space_by_index ( mesh_index, 1 ) );
        }

        this->mesh_list_ = new_mesh_list;
        this->mesh_change_times_.data ( ) = new_change_list.data ( );

        std::vector<int> new_mesh_indices = indicator_mesh_indices;
        for ( int t = 0; t < new_mesh_indices.size ( ); ++t )
        {
            int old_index = indicator_mesh_indices[t];
            int pot_new_index = this->mesh_index ( this->time_mesh_->time ( t ) );
            int pos = -1;
            if ( old2new_mesh_indices[old_index].find ( pot_new_index, &pos ) )
            {
                new_mesh_indices[t] = pot_new_index;
            }
            else
            {
                new_mesh_indices[t] = old2new_mesh_indices[old_index][0];
            }
        }

        LOG_DEBUG ( 2, " old indicator 2 mesh index: " << string_from_range ( indicator_mesh_indices.begin ( ), indicator_mesh_indices.end ( ) ) );

        indicator_mesh_indices = new_mesh_indices;

        LOG_DEBUG ( 2, " new indicator 2 mesh index: " << string_from_range ( indicator_mesh_indices.begin ( ), indicator_mesh_indices.end ( ) ) );

        /* write mesh list */
    }

    template<class LAD, IMPL MESH_IMPL, int DIM>
    void DynamicMeshHandler<LAD, MESH_IMPL, DIM>::refine ( int mesh_index, std::vector<int>& markers )
    {

        MeshPtr active_mesh = this->get_mesh_by_index ( mesh_index );
        if ( MESH_IMPL == mesh::IMPL_P4EST )
        {
            boost::intrusive_ptr<MeshPXest> mesh_pXest = boost::static_pointer_cast<MeshPXest> ( active_mesh );
            mesh_pXest->set_patch_mode ( true );
            mesh_pXest->set_connection_mode ( 2 );
        }

        SharedVertexTable shared_verts;

        active_mesh = active_mesh->refine ( markers );
        if ( MESH_IMPL == IMPL_P4EST )
        {
            this->mesh_list_[mesh_index] = compute_ghost_cells ( *active_mesh, this->comm_, shared_verts, mesh::IMPL_P4EST, 2 );
        }
        else if ( MESH_IMPL == IMPL_DBVIEW )
        {
            this->mesh_list_[mesh_index] = compute_ghost_cells ( *active_mesh, this->comm_, shared_verts, mesh::IMPL_DBVIEW, 1 );
        }
        active_mesh = this->get_mesh_by_index ( mesh_index );
        int num_global_cells = active_mesh->num_global_cells ( this->comm_ );

        LOG_DEBUG ( 1, "[" << rank_ << "] "
                    << "  Mesh " << mesh_index
                    << "  #local cells: " << active_mesh->num_local_cells ( )
                    << ", #ghost cells: " << active_mesh->num_ghost_cells ( )
                    << ", #sum: " << active_mesh->num_entities ( DIM )
                    << ", #global: " << num_global_cells );

    }

    template<class LAD, IMPL MESH_IMPL, int DIM>
    void DynamicMeshHandler<LAD, MESH_IMPL, DIM>::load_mesh_from_file ( int adapt_counter, int num_mesh, std::string& prefix )
    {
        this->mesh_list_.clear ( );
        this->mesh_list_.resize ( num_mesh );

        for ( int l = 0; l < num_mesh; ++l )
        {
            std::stringstream pre;
            pre << prefix << "." << adapt_counter << "." << l << ".h5";
            std::string filename = pre.str ( );

            this->mesh_list_[l] = load_mesh ( filename, this->comm_, MESH_IMPL );
        }
        this->problem_->set_active_mesh ( this->mesh_list_[0] );
    }

    template<class LAD, IMPL MESH_IMPL, int DIM>
    void DynamicMeshHandler<LAD, MESH_IMPL, DIM>::visualize_mesh ( int adapt_counter, std::string& prefix )
    {
        for ( int m = 0; m<this->num_mesh ( ); ++m )
        {
            // visualize adapt markers
            /*
            std::vector<double> attr_val(num_cells, 0.);
            for (int c=0; c<num_cells; ++c)
            {
            attr_val[c] = adapt_markers[c];
            }
            AttributePtr attr_est ( new DoubleAttribute ( attr_val ) );
            active_mesh->add_attribute ( "marker" , DIM, attr_est );

            std::stringstream a_input;
            std::stringstream a_pre;
            a_pre << this->root_ << "/mesh/adapt_mesh." << this->adapt_counter_ << "." << m;

            if(this->num_partitions_ > 1)
            a_input << ".pvtu";
            else 
            a_input << ".vtu";

            std::string a_filename = a_pre.str() + a_input.str();

            PVtkWriter a_writer ( this->comm_ );
            std::ostringstream a_name;
            a_name << this->root_ << "/mesh/adapt_mesh." << this->adapt_counter_ << "." << m << ".pvtu";
            std::string a_output_file = a_name.str ( );
            a_writer.add_all_attributes ( *active_mesh, true );
            a_writer.write ( a_output_file.c_str ( ), *active_mesh );
             */

            PVtkWriter writer ( this->comm_ );
            std::ostringstream name;
            name << prefix << "." << adapt_counter << "." << m << ".pvtu";
            std::string output_file = name.str ( );
            writer.add_all_attributes ( *this->get_mesh_by_index ( m ), true );
            writer.write ( output_file.c_str ( ), *this->get_mesh_by_index ( m ) );
        }
    }

    template<class LAD, IMPL MESH_IMPL, int DIM>
    void DynamicMeshHandler<LAD, MESH_IMPL, DIM>::save_mesh_to_file ( int adapt_counter, std::string& prefix )
    {
        for ( int l = 0; l<this->mesh_list_.size ( ); ++l )
        {
            std::stringstream pre;
            pre << prefix << "." << adapt_counter << "." << l << ".h5";
            std::string filename = pre.str ( );

            save_mesh ( filename, this->mesh_list_[l], this->comm_ );
        }
    }

    template class DynamicMeshHandler<LADescriptorCoupledD, IMPL_P4EST, 1>;
    template class DynamicMeshHandler<LADescriptorCoupledD, IMPL_P4EST, 2>;
    template class DynamicMeshHandler<LADescriptorCoupledD, IMPL_P4EST, 3>;
    template class DynamicMeshHandler<LADescriptorCoupledD, IMPL_DBVIEW, 1>;
    template class DynamicMeshHandler<LADescriptorCoupledD, IMPL_DBVIEW, 2>;
    template class DynamicMeshHandler<LADescriptorCoupledD, IMPL_DBVIEW, 3>;

#ifdef WITH_HYPRE
    template class DynamicMeshHandler<LADescriptorHypreD, IMPL_DBVIEW, 1>;
    template class DynamicMeshHandler<LADescriptorHypreD, IMPL_DBVIEW, 2>;
    template class DynamicMeshHandler<LADescriptorHypreD, IMPL_DBVIEW, 3>;

    template class DynamicMeshHandler<LADescriptorHypreD, IMPL_P4EST, 1>;
    template class DynamicMeshHandler<LADescriptorHypreD, IMPL_P4EST, 2>;
    template class DynamicMeshHandler<LADescriptorHypreD, IMPL_P4EST, 3>;
#endif

    //////////////////////////////////////////////////////////////////
    //////////// DynamicMeshProblem //////////////////////////////////

    template<class LAD, IMPL MESH_IMPL, int DIM>
    void DynamicMeshProblem<LAD, MESH_IMPL, DIM>::init_dmh ( std::string& adapt_type )
    {
        // pass initial time mesh to handler
        this->dmh_->set_time_mesh ( &this->t_mesh_ );

        // set type mesh change process
        this->dmh_->set_adapt_type ( adapt_type );

        // pass vectors to handler which have to be updated for each mesh change
        this->set_update_vectors ( );

        // set initial mesh change points and pass to handler
        this->init_mesh_change_list ( );
    }

    template<class LAD, IMPL MESH_IMPL, int DIM>
    void DynamicMeshProblem<LAD, MESH_IMPL, DIM>::interpolate_vector ( const typename LAD::VectorType& in_vec, int in_index,
                                                                       typename LAD::VectorType& out_vec, int out_index,
                                                                       int mode )
    {
        if ( in_index == -1 )
        {
            return;
        }
        if ( in_index == out_index )
        {
            out_vec.CopyFrom ( in_vec );
            return;
        }
        this->dmh_->interpolate_vector ( in_vec, in_index, out_vec, out_index );
    }

    template<class LAD, IMPL MESH_IMPL, int DIM>
    void DynamicMeshProblem<LAD, MESH_IMPL, DIM>::adapt_mesh_change_list ( std::string adapt_type, int adapt_counter,
                                                                           int min_steps_for_mesh, int max_mesh_number,
                                                                           int rate, typename LAD::DataType tol,
                                                                           const std::vector< std::vector< typename LAD::DataType> >& space_indicator,
                                                                           std::vector<int>& indicator_mesh_indices )
    {
        if ( adapt_type == "None" || adapt_type == "Fixed" )
        {
            return;
        }

        int rank;
        MPI_Comm_rank ( this->dmh_->get_mpi_comm ( ), &rank );

        int num_steps = this->t_mesh_.num_intervals ( adapt_counter );
        int num_mesh = this->dmh_->num_mesh ( );

        // Compute differences of error estimators for each time step
        std::vector<double> global_diff ( num_steps - 1, 0. );
        std::vector<double> local_diff ( num_steps - 1, 0. );
        std::vector<double> local_norm ( num_steps - 1, 0. );
        std::vector<double> global_norm ( num_steps - 1, 0. );

        // Loop over all meshes
        for ( int m = 0; m < num_mesh; ++m )
        {
            int first_t = this->dmh_->first_step_for_mesh ( m );
            int last_t = this->dmh_->last_step_for_mesh ( m );

            int num_cells = this->dmh_->num_cells_by_step ( first_t );

            // Loop over all time step belongin to current mesh, i.e. we ignore differences between time steps belonging to different meshes
            for ( int t = first_t; t < last_t; ++t )
            {
                // Loop over all cells in current mesh
                for ( int c = 0; c < num_cells; ++c )
                {
                    // get space indicator for current cell c and timestep t
                    double old_val = space_indicator[t][c];
                    double new_val = space_indicator[t + 1][c];

                    local_diff[t] += ( old_val - new_val ) * ( old_val - new_val );
                    local_norm[t] += old_val * old_val;
                }
            }
        }

        // Allreduce
        MPI_Allreduce ( &local_diff[0], &global_diff[0], num_steps - 1, MPI_DOUBLE, MPI_SUM, this->dmh_->get_mpi_comm ( ) );
        MPI_Allreduce ( &local_norm[0], &global_norm[0], num_steps - 1, MPI_DOUBLE, MPI_SUM, this->dmh_->get_mpi_comm ( ) );

        std::vector<double> rel_diff ( num_steps - 1, 0. );

        for ( int l = 0; l < global_diff.size ( ); ++l )
        {
            global_diff[l] = std::sqrt ( global_diff[l] );
            global_norm[l] = std::sqrt ( global_norm[l] );

            rel_diff[l] = global_diff[l] / global_norm[l];
        }

        std::vector<int> change_steps;

        // sort differences in ascending order
        std::vector<int> sort_ind ( global_diff.size ( ), 0 );
        for ( int i = 0; i < rel_diff.size ( ); ++i )
        {
            sort_ind[i] = i;
        }
        sortingPermutation ( rel_diff, sort_ind );

        SortedArray<int> mesh_change_steps;
        mesh_change_steps.data ( ) = this->dmh_->get_mesh_change_steps ( );

        LOG_DEBUG ( 1, "Old mesh change steps " << string_from_range ( mesh_change_steps.begin ( ), mesh_change_steps.end ( ) ) );

        // insert fixed number of mesh change points 
        if ( adapt_type == "FixedFraction" )
        {
            // find #rate largest differences
            int added_steps = 0;
            for ( int t = sort_ind.size ( ) - 1; t >= 0; --t )
            {
                if ( added_steps >= rate )
                {
                    break;
                }
                int trial_step = sort_ind[t];

                // check if maximum numbre of mesh changes is reached 
                if ( change_steps.size ( ) + num_mesh < max_mesh_number )
                {
                    // check if trial step is not first mesh
                    if ( trial_step > 0 )
                    {
                        int pos = -1;
                        // check if trial step is already contained in list of mesh change steps
                        if ( !mesh_change_steps.find ( trial_step, &pos ) )
                        {
                            mesh_change_steps.insert ( trial_step );
                            mesh_change_steps.find ( trial_step, &pos );
                            int prev_step = -1;
                            int next_step = -1;
                            if ( pos > 0 )
                            {
                                prev_step = mesh_change_steps.data ( ).at ( pos - 1 );
                            }
                            if ( pos < mesh_change_steps.size ( ) - 1 )
                            {
                                next_step = mesh_change_steps.data ( ).at ( pos + 1 );
                            }

                            // check if distance to previuous change is large enough 
                            if ( prev_step >= 0 && ( trial_step - prev_step ) <= min_steps_for_mesh )
                            {
                                mesh_change_steps.erase ( pos );
                                continue;
                            }

                            // check if distance to next change step is large enough
                            if ( next_step >= 0 && ( next_step - trial_step ) <= min_steps_for_mesh )
                            {
                                mesh_change_steps.erase ( pos );
                                continue;
                            }

                            // accept trial step
                            change_steps.push_back ( trial_step );
                            added_steps++;
                        }
                    }
                }
            }
        }
            // insert mesh change points according to given threshold 
        else if ( adapt_type == "FixedError" )
        {
            // loop over differences is descending order
            for ( int t = sort_ind.size ( ) - 1; t >= 0; --t )
            {
                int trial_step = sort_ind[t];

                // check if change of trial_step is large enough
                if ( rel_diff[trial_step] >= tol )
                {
                    // check if maximum number of mesh change steps is reached
                    if ( change_steps.size ( ) + num_mesh < max_mesh_number )
                    {
                        // check if trial step is not first step
                        if ( trial_step > 0 )
                        {
                            int pos = -1;
                            // check if trial step is already contained in list of change steps
                            if ( !mesh_change_steps.find ( trial_step, &pos ) )
                            {
                                mesh_change_steps.insert ( trial_step );
                                mesh_change_steps.find ( trial_step, &pos );
                                int prev_step = -1;
                                int next_step = -1;
                                if ( pos > 0 )
                                {
                                    prev_step = mesh_change_steps.data ( ).at ( pos - 1 );
                                }
                                if ( pos < mesh_change_steps.size ( ) - 1 )
                                {
                                    next_step = mesh_change_steps.data ( ).at ( pos + 1 );
                                }

                                // check if distance to previuous change is large enough 
                                if ( prev_step >= 0 && ( trial_step - prev_step ) <= min_steps_for_mesh )
                                {
                                    mesh_change_steps.erase ( pos );
                                    continue;
                                }

                                // check if distance to next change step is large enough
                                if ( next_step >= 0 && ( next_step - trial_step ) <= min_steps_for_mesh )
                                {
                                    mesh_change_steps.erase ( pos );
                                    continue;
                                }

                                // accept trial step
                                change_steps.push_back ( trial_step );
                            }
                        }
                    }
                }
            }
        }

        LOG_DEBUG ( 1, "New mesh change steps " << string_from_range ( mesh_change_steps.begin ( ), mesh_change_steps.end ( ) ) );

        // adapt dynamic mesh handler
        std::vector<double> add_times;
        for ( int t = 0; t < change_steps.size ( ); ++t )
        {
            add_times.push_back ( this->t_mesh_.time ( change_steps[t], adapt_counter ) );
        }
        this->dmh_->add_mesh_change_times ( add_times, indicator_mesh_indices );
    }

    template class DynamicMeshProblem<LADescriptorCoupledD, IMPL_P4EST, 1>;
    template class DynamicMeshProblem<LADescriptorCoupledD, IMPL_P4EST, 2>;
    template class DynamicMeshProblem<LADescriptorCoupledD, IMPL_P4EST, 3>;
    template class DynamicMeshProblem<LADescriptorCoupledD, IMPL_DBVIEW, 1>;
    template class DynamicMeshProblem<LADescriptorCoupledD, IMPL_DBVIEW, 2>;
    template class DynamicMeshProblem<LADescriptorCoupledD, IMPL_DBVIEW, 3>;

#ifdef WITH_HYPRE
    template class DynamicMeshProblem<LADescriptorHypreD, IMPL_DBVIEW, 1>;
    template class DynamicMeshProblem<LADescriptorHypreD, IMPL_DBVIEW, 2>;
    template class DynamicMeshProblem<LADescriptorHypreD, IMPL_DBVIEW, 3>;

    template class DynamicMeshProblem<LADescriptorHypreD, IMPL_P4EST, 1>;
    template class DynamicMeshProblem<LADescriptorHypreD, IMPL_P4EST, 2>;
    template class DynamicMeshProblem<LADescriptorHypreD, IMPL_P4EST, 3>;
#endif
}
