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

/// \author Simon Gawlok

#include "preconditioner_vanka.h"

namespace hiflow
{
    namespace la
    {

        template<class LAD>
        PreconditionerVanka<LAD>::PreconditionerVanka ( )
        : PreconditionerBlockJacobi<LAD>( )
        {
        }

        template<class LAD>
        PreconditionerVanka<LAD>::~PreconditionerVanka ( )
        {
            this->Clear ( );
        }

        template<class LAD>
        void PreconditionerVanka<LAD>::InitParameter
        (
          const hiflow::VectorSpace<DataType>& space,
          const DataType damping_param,
          const int num_iter,
          const bool use_block_of_cells,
          const bool use_preconditioner
          )
        {
            this->space_ = &space;
            this->damping_param_ = damping_param;
            this->maxits_ = num_iter;
            this->use_block_of_cells_ = use_block_of_cells;
            this->use_preconditioner_ = use_preconditioner;
            if ( this->print_level_ > 2 )
            {
                LOG_INFO ( "Damping parameter", this->damping_param_ );
                LOG_INFO ( "Number of iterations", this->maxits_ );
            }
        }

        template<class LAD>
        void PreconditionerVanka<LAD>::InitIluppPrecond
        (
          int prepro_type,
          int precond_no,
          int max_levels,
          DataType mem_factor,
          DataType threshold,
          DataType min_pivot
          )
        {
            this->precond_.InitParameter ( prepro_type,
                                           precond_no,
                                           max_levels,
                                           mem_factor,
                                           threshold,
                                           min_pivot );
            this->use_preconditioner_ = true;
        }

        template<class LAD>
        void PreconditionerVanka<LAD>::Factorize ( )
        {
            if ( !this->use_block_of_cells_ )
            {
                for ( int i = 0; i < local_mat_diag_.size ( ); ++i )
                {
                    local_mat_diag_[i]->Factorize ( );
                }
            }
        }

        template<class LAD>
        void PreconditionerVanka<LAD>::SetupOperator ( OperatorType& op )
        {
            this->op_ = &op;
            this->SetModifiedOperator ( true );

            if ( this->use_preconditioner_ )
            {
                this->precond_.SetupOperator ( *( this->op_ ) );
                this->precond_.SetModifiedOperator ( true );
            }
        }

        template<class LAD>
        void PreconditionerVanka<LAD>::Build ( )
        {
            assert ( this->op_ != NULL );
            if ( this->print_level_ > 2 )
            {
                LOG_INFO ( "Build Solver", 1 );
            }

            this->CreateLocalMatrices ( );
            this->Factorize ( );

            this->SetState ( true );
            this->SetModifiedOperator ( false );

            if ( this->use_preconditioner_ )
            {
                if ( !this->precond_.GetReuse ( ) || !this->precond_.GetState ( ) )
                {
                    this->precond_.Build ( );
                    this->precond_.SetState ( true );
                    this->precond_.SetModifiedOperator ( false );
                }
            }
        }

        template<class LAD>
        LinearSolverState PreconditionerVanka<LAD>::ApplyPreconditioner
        (
          const VectorType& b,
          VectorType* x )
        {

            if ( !this->GetState ( ) )
            {
                this->Build ( );
            }

            if ( this->use_preconditioner_ )
            {
                this->precond_.ApplyPreconditioner ( b, x );
            }

            // locally needed objects
            std::vector<DataType> b_loc;
            std::vector<DataType> x_loc;
            std::vector<DataType> x_temp_loc;
            std::vector<DataType> res_loc;
            std::vector<DataType> res_loc_2;

            for ( int k = 0, k_e = this->maxits_; k != k_e; ++k )
            {
                x->Update ( );

                // Loop over all "cells"
                for ( size_t i = 0, i_e = sorted_dofs_diag_.size ( ); i != i_e; ++i )
                {

                    // COMPUTE CURRENT RESIDUAL

                    res_loc.resize
                            ( sorted_dofs_diag_[i].size ( ),
                              0. );
                    res_loc_2.resize
                            ( sorted_dofs_diag_[i].size ( ),
                              0. );
                    x_temp_loc.resize
                            ( sorted_dofs_diag_[i].size ( ),
                              0. );

                    this->op_->VectorMult_submatrix_vanka
                            ( vec2ptr ( sorted_dofs_diag_[i] ),
                              sorted_dofs_diag_[i].size ( ),
                              *x,
                              vec2ptr ( res_loc ) );

                    x->GetValues ( vec2ptr ( sorted_dofs_diag_[i] ),
                                   sorted_dofs_diag_[i].size ( ),
                                   vec2ptr ( x_temp_loc ) );

                    if ( !this->use_block_of_cells_ )
                    {
                        local_mat_diag_[i]->VectorMult ( x_temp_loc, res_loc_2 );
                    }
                    else
                    {
                        // Resize local matrix to correct size
                        local_mat_diag_block_mode_.Resize
                                ( sorted_dofs_diag_[i].size ( ),
                                  sorted_dofs_diag_[i].size ( )
                                  );

                        // Get local matrix from global operator
                        this->op_->GetValues
                                (
                                  vec2ptr ( sorted_dofs_diag_[i] ),
                                  sorted_dofs_diag_[i].size ( ),
                                  vec2ptr ( sorted_dofs_diag_[i] ),
                                  sorted_dofs_diag_[i].size ( ),
                                  &( ( local_mat_diag_block_mode_ ) ( 0, 0 ) )
                                  );

                        // Set block-size for LU decomposition
                        local_mat_diag_block_mode_.set_blocksize ( 16 );

                        local_mat_diag_block_mode_.VectorMult ( x_temp_loc, res_loc_2 );
                    }

                    // Get local rhs
                    b_loc.resize
                            ( sorted_dofs_diag_[i].size ( ),
                              0. );
                    b.GetValues
                            ( vec2ptr ( sorted_dofs_diag_[i] ),
                              sorted_dofs_diag_[i].size ( ),
                              vec2ptr ( b_loc ) );

                    for ( size_t j = 0, j_e = b_loc.size ( ); j != j_e; ++j )
                    {
                        b_loc[j] = b_loc[j] - res_loc[j] + res_loc_2[j];
                    }

                    // Solve local problem
                    x_loc.resize
                            ( sorted_dofs_diag_[i].size ( ),
                              0. );
                    if ( !this->use_block_of_cells_ )
                    {
                        local_mat_diag_[i]->ForwardBackward
                                ( b_loc, x_loc );
                    }
                    else
                    {
                        local_mat_diag_block_mode_.Solve ( b_loc, x_loc );
                    }

                    // add local contribution to solution vector
                    for ( size_t j = 0, j_e = x_loc.size ( ); j != j_e; ++j )
                    {
                        x_temp_loc[j] =
                                ( 1. - this->damping_param_ )
                                * x_temp_loc[j]
                                + this->damping_param_
                                * x_loc[j];
                    }

                    x->SetValues ( vec2ptr ( sorted_dofs_diag_[i] ),
                                   sorted_dofs_diag_[i].size ( ),
                                   vec2ptr ( x_temp_loc ) );
                }
                // Loop backward over all "cells"
                for ( int i = sorted_dofs_diag_.size ( ); i--; )
                {

                    // COMPUTE CURRENT RESIDUAL

                    res_loc.resize
                            ( sorted_dofs_diag_[i].size ( ),
                              0. );
                    res_loc_2.resize
                            ( sorted_dofs_diag_[i].size ( ),
                              0. );
                    x_temp_loc.resize
                            ( sorted_dofs_diag_[i].size ( ),
                              0. );

                    this->op_->VectorMult_submatrix_vanka
                            ( vec2ptr ( sorted_dofs_diag_[i] ),
                              sorted_dofs_diag_[i].size ( ),
                              *x,
                              vec2ptr ( res_loc ) );

                    x->GetValues ( vec2ptr ( sorted_dofs_diag_[i] ),
                                   sorted_dofs_diag_[i].size ( ),
                                   vec2ptr ( x_temp_loc ) );

                    if ( !this->use_block_of_cells_ )
                    {
                        local_mat_diag_[i]->VectorMult ( x_temp_loc, res_loc_2 );
                    }
                    else
                    {
                        // Resize local matrix to correct size
                        local_mat_diag_block_mode_.Resize
                                ( sorted_dofs_diag_[i].size ( ),
                                  sorted_dofs_diag_[i].size ( )
                                  );

                        // Get local matrix from global operator
                        this->op_->GetValues
                                (
                                  vec2ptr ( sorted_dofs_diag_[i] ),
                                  sorted_dofs_diag_[i].size ( ),
                                  vec2ptr ( sorted_dofs_diag_[i] ),
                                  sorted_dofs_diag_[i].size ( ),
                                  &( ( local_mat_diag_block_mode_ ) ( 0, 0 ) )
                                  );

                        // Set block-size for LU decomposition
                        local_mat_diag_block_mode_.set_blocksize ( 16 );

                        local_mat_diag_block_mode_.VectorMult ( x_temp_loc, res_loc_2 );
                    }

                    // Get local rhs
                    b_loc.resize
                            ( sorted_dofs_diag_[i].size ( ),
                              0. );
                    b.GetValues
                            ( vec2ptr ( sorted_dofs_diag_[i] ),
                              sorted_dofs_diag_[i].size ( ),
                              vec2ptr ( b_loc ) );

                    for ( size_t j = 0, j_e = b_loc.size ( ); j != j_e; ++j )
                    {
                        b_loc[j] = b_loc[j] - res_loc[j] + res_loc_2[j];
                    }

                    // Solve local problem
                    x_loc.resize
                            ( sorted_dofs_diag_[i].size ( ),
                              0. );
                    if ( !this->use_block_of_cells_ )
                    {
                        local_mat_diag_[i]->ForwardBackward
                                ( b_loc, x_loc );
                    }
                    else
                    {
                        local_mat_diag_block_mode_.Solve ( b_loc, x_loc );
                    }

                    // add local contribution to solution vector
                    for ( size_t j = 0, j_e = x_loc.size ( ); j != j_e; ++j )
                    {
                        x_temp_loc[j] =
                                ( 1. - this->damping_param_ )
                                * x_temp_loc[j]
                                + this->damping_param_
                                * x_loc[j];
                    }

                    x->SetValues ( vec2ptr ( sorted_dofs_diag_[i] ),
                                   sorted_dofs_diag_[i].size ( ),
                                   vec2ptr ( x_temp_loc ) );
                }
            }
            x->Update ( );

            return kSolverSuccess;
        }

        template<class LAD>
        void PreconditionerVanka<LAD>::Clear ( )
        {

            if ( !this->use_block_of_cells_ )
            {
                for ( size_t i = 0, i_e = sorted_dofs_diag_.size ( ); i != i_e; ++i )
                {
                    local_mat_diag_[i]->Clear ( );
                    delete local_mat_diag_[i];
                    sorted_dofs_diag_[i].clear ( );
                }
                local_mat_diag_.clear ( );
            }
            else
            {
                this->local_mat_diag_block_mode_.Clear ( );
            }

            sorted_dofs_diag_.clear ( );
            if ( this->use_preconditioner_ )
            {
                this->precond_.Clear ( );
            }
            Preconditioner<LAD>::Clear ( );
        }

        template<class LAD>
        void PreconditionerVanka<LAD>::CreateLocalMatrices ( )
        {
            // Get pointer to finite element mesh
            hiflow::mesh::ConstMeshPtr mesh = &( space_->mesh ( ) );

            // Get number of elements in (local) mesh
            const int num_elements = mesh->num_entities ( mesh->tdim ( ) );

            if ( !this->use_block_of_cells_ )
            {
                // prepare local matrices
                local_mat_diag_.resize ( num_elements, new SeqDenseMatrix<DataType> );
                for ( int i = 0, i_e = local_mat_diag_.size ( ); i != i_e; ++i )
                {
                    local_mat_diag_[i] = new SeqDenseMatrix<DataType>;
                    local_mat_diag_[i]->Resize ( 0, 0 );
                }
            }

            // prepare indices vectors
            sorted_dofs_diag_.resize ( num_elements );

            // Get local matrices from global matrix
            for ( hiflow::mesh::EntityIterator it = mesh->begin ( mesh->tdim ( ) ), e_it = mesh->end ( mesh->tdim ( ) );
                  it != e_it;
                  ++it )
            {
                // get global dof indices on current cell
                for ( int var = 0, var_e = space_->get_nb_var ( ); var != var_e; ++var )
                {
                    std::vector<int> var_indices;
                    space_->dof ( ).get_dofs_on_cell ( var, it->index ( ), var_indices );
                    sorted_dofs_diag_[it->index ( )].insert
                            ( sorted_dofs_diag_[it->index ( )].end ( ),
                              var_indices.begin ( ),
                              var_indices.end ( ) );
                }

                // get dofs of all adjacent cells in block-mode
                if ( this->use_block_of_cells_ )
                {
                    for ( hiflow::mesh::IncidentEntityIterator it_inc = it->begin_incident ( mesh->tdim ( ) ),
                          e_it_inc = it->end_incident ( mesh->tdim ( ) );
                          it_inc != e_it_inc;
                          ++it_inc )
                    {
                        // get global dof indices on current incident cell
                        for ( int var = 0, var_e = space_->get_nb_var ( ); var != var_e; ++var )
                        {
                            std::vector<int> var_indices;
                            space_->dof ( ).get_dofs_on_cell ( var, it_inc->index ( ), var_indices );
                            sorted_dofs_diag_[it->index ( )].insert
                                    ( sorted_dofs_diag_[it->index ( )].end ( ),
                                      var_indices.begin ( ),
                                      var_indices.end ( ) );
                        }
                    }
                }

                // sort dofs and make them unique
                std::sort ( sorted_dofs_diag_[it->index ( )].begin ( ),
                            sorted_dofs_diag_[it->index ( )].end ( ) );
                sorted_dofs_diag_[it->index ( )].erase
                        (
                          unique
                          (
                            sorted_dofs_diag_[it->index ( )].begin ( ),
                            sorted_dofs_diag_[it->index ( )].end ( )
                            ),
                          sorted_dofs_diag_[it->index ( )].end ( )
                          );

                // reduce dof indices to only subdomain indices
                bool found = true;
                while ( found )
                {
                    found = false;
                    for ( int i = 0, i_e = sorted_dofs_diag_[it->index ( )].size ( );
                          i != i_e && !found;
                          ++i )
                    {
                        if (
                             !this->space_->dof ( ).is_dof_on_sd
                             (
                               sorted_dofs_diag_[it->index ( )][i]
                               )
                             )
                        {
                            sorted_dofs_diag_[it->index ( )].erase
                                    (
                                      sorted_dofs_diag_[it->index ( )].begin ( ) + i
                                      );
                            found = true;
                        }
                    }
                }
            }

            // remove cells with empty dof list
            bool found = true;
            while ( found )
            {
                found = false;
                for ( int i = 0, i_e = sorted_dofs_diag_.size ( );
                      i != i_e && !found;
                      ++i )
                {
                    if ( sorted_dofs_diag_[i].size ( ) == 0
                         )
                    {
                        if ( !this->use_block_of_cells_ )
                        {
                            local_mat_diag_[i]->Clear ( );
                            delete local_mat_diag_[i];
                            local_mat_diag_.erase
                                    (
                                      local_mat_diag_.begin ( ) + i
                                      );
                        }
                        sorted_dofs_diag_.erase
                                (
                                  sorted_dofs_diag_.begin ( ) + i
                                  );
                        found = true;
                    }
                }
            }

            // Get local matrices
            if ( !this->use_block_of_cells_ )
            {
                for ( int i = 0, i_e = local_mat_diag_.size ( ); i != i_e; ++i )
                {
                    // Resize local matrix to correct size
                    local_mat_diag_[i]->Resize
                            ( sorted_dofs_diag_[i].size ( ),
                              sorted_dofs_diag_[i].size ( )
                              );

                    // Get local matrix from global operator
                    this->op_->GetValues
                            (
                              vec2ptr ( sorted_dofs_diag_[i] ),
                              sorted_dofs_diag_[i].size ( ),
                              vec2ptr ( sorted_dofs_diag_[i] ),
                              sorted_dofs_diag_[i].size ( ),
                              &( ( *local_mat_diag_[i] )( 0, 0 ) )
                              );

                    // Set block-size for LU decomposition
                    local_mat_diag_[i]->set_blocksize ( 16 );
                }
            }
        }

        /// template instantiation
        template class PreconditionerVanka<hiflow::la::LADescriptorCoupledD>;
        //template class PreconditionerVanka<hiflow::la::LADescriptorCoupledS>;

    } // namespace la
} // namespace hiflow
