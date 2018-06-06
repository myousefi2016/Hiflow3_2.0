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

/// @author Philipp Gerstner

#include <nonlinear/nonlinear_problem.h>

#include "tri_block.h"

namespace hiflow
{
    namespace la
    {

        template<class LAD>
        TriBlock<LAD>::TriBlock ( )
        : SchurComplement<LAD> ( )
        {
            this->precond_ = NULL;
            this->standard_precond_ = false;

            this->solver_E_ = NULL;
            this->scaling_E_ = 0.0;

            this->block_op_.resize ( 4, 0 );
            this->E_is_initialized_ = false;
            this->E_is_set_ = false;
            this->E_modified_ = false;
            this->E_passed2solver_ = false;
        }

        template<class LAD>
        TriBlock<LAD>::~TriBlock ( )
        {
            this->Clear ( );
        }

        template<class LAD>
        void TriBlock<LAD>::Clear ( )
        {
            SchurComplement<LAD>::Clear ( );

            this->E_.Clear ( );
            this->E_is_initialized_ = false;
            this->E_is_set_ = false;
            this->diag_D_.clear ( );

            if ( this->solver_E_ != NULL )
            {
                this->solver_E_->Clear ( );
                this->solver_E_ = NULL;
            }

            this->scaling_E_ = 0.0;
            this->scaling_S_ = 0.0;
            this->block_op_.resize ( 4, 0 );
            this->standard_precond_ = false;
        }

        template<class LAD>
        void TriBlock<LAD>::Init
        ( const hiflow::VectorSpace<ValueType>& space,
          const std::vector<int>& block_one_variables,
          const std::vector<int>& block_two_variables
          )
        {
            SchurComplement<LAD>::Init ( space, block_one_variables, block_two_variables );

            E_.Init ( this->comm_, this->la_c_two_, this->la_c_two_ );
        }

        template<class LAD>
        void TriBlock<LAD>::SetupOperator
        ( OperatorType& op, const int block_id )
        {
            if ( this->use_press_conv_diff_ && block_id == 3 )
            {
                std::cout << "TriBlock: Don't SetupOperator D if you want to use ApproximateSchur. " << std::endl;
                exit ( -1 );
            }

            //*****************************************************************
            // Operator setup means copying entries from system matrix op to
            // submatrices A, B, C and D
            //*****************************************************************
            this->block_op_[block_id] = 1;

            // Operator / Block A
            if ( block_id == 0 )
            {
                this->A_.Zeros ( );

                int row_index_block = 0;
                for ( typename std::map<int, SortedArray<int> >::const_iterator it = this->couplings_A_.begin ( ),
                      e_it = this->couplings_A_.end ( );
                      it != e_it;
                      ++it )
                {
                    // Get values from system matrix op
                    int row_index_system = it->first;
                    std::vector<ValueType> vals ( it->second.size ( ), 0. );
                    op.GetValues ( &row_index_system, 1, vec2ptr ( it->second.data ( ) ), it->second.size ( ), vec2ptr ( vals ) );

                    // Set values in glock matrix A
                    int block_index_system = row_index_block + this->offsets_block_one_[this->my_rank_];
                    this->A_.SetValues ( &block_index_system, 1, vec2ptr ( this->sparsity_A_[row_index_block] ), it->second.size ( ), vec2ptr ( vals ) );
                    ++row_index_block;
                }
                this->A_modified_ = true;
                this->A_passed2solver_ = false;
            }

            // sub matrix B
            if ( block_id == 1 )
            {
                this->B_.Zeros ( );
                int row_index_block = 0;
                for ( typename std::map<int, SortedArray<int> >::const_iterator it = this->couplings_B_.begin ( ),
                      e_it = this->couplings_B_.end ( );
                      it != e_it;
                      ++it )
                {
                    // Get values from system matrix op
                    int row_index_system = it->first;
                    std::vector<ValueType> vals ( it->second.size ( ), 0. );
                    op.GetValues ( &row_index_system, 1, vec2ptr ( it->second.data ( ) ), it->second.size ( ), vec2ptr ( vals ) );

                    // Set values in glock matrix B
                    int block_index_system = row_index_block + this->offsets_block_one_[this->my_rank_];
                    this->B_.SetValues ( &block_index_system, 1, vec2ptr ( this->sparsity_B_[row_index_block] ), it->second.size ( ), vec2ptr ( vals ) );
                    ++row_index_block;
                }
            }

            // sub matrix C
            if ( block_id == 2 )
            {
                this->C_.Zeros ( );
                int row_index_block = 0;
                for ( typename std::map<int, SortedArray<int> >::const_iterator it = this->couplings_C_.begin ( ),
                      e_it = this->couplings_C_.end ( );
                      it != e_it;
                      ++it )
                {
                    // Get values from system matrix op
                    int row_index_system = it->first;
                    std::vector<ValueType> vals ( it->second.size ( ), 0. );
                    op.GetValues ( &row_index_system, 1, vec2ptr ( it->second.data ( ) ), it->second.size ( ), vec2ptr ( vals ) );

                    // Set values in block matrix C
                    int block_index_system = row_index_block + this->offsets_block_two_[this->my_rank_];
                    this->C_.SetValues ( &block_index_system, 1, vec2ptr ( this->sparsity_C_[row_index_block] ), it->second.size ( ), vec2ptr ( vals ) );
                    ++row_index_block;
                }
            }

            // sub matrix D
            if ( block_id == 3 )
            {
                this->D_.Zeros ( );
                int row_index_block = 0;
                for ( typename std::map<int, SortedArray<int> >::const_iterator it = this->couplings_D_.begin ( ),
                      e_it = this->couplings_D_.end ( );
                      it != e_it;
                      ++it )
                {
                    // Get values from system matrix op
                    int row_index_system = it->first;
                    std::vector<ValueType> vals ( it->second.size ( ), 0. );
                    op.GetValues ( &row_index_system, 1, vec2ptr ( it->second.data ( ) ), it->second.size ( ), vec2ptr ( vals ) );

                    // Set values in block matrix D
                    int block_index_system = row_index_block + this->offsets_block_two_[this->my_rank_];
                    this->D_.SetValues ( &block_index_system, 1, vec2ptr ( this->sparsity_D_[row_index_block] ), it->second.size ( ), vec2ptr ( vals ) );
                    ++row_index_block;
                }
                this->D_modified_ = true;
                this->D_passed2solver_ = false;
            }

            //*****************************************************************
            // Initialize diag_A_
            //*****************************************************************
            if ( block_id == 0 )
            {
                this->diag_A_.clear ( );
                this->diag_A_.resize ( this->sparsity_A_.size ( ), 0. );

                const int offset_A = this->offsets_block_one_[this->my_rank_];
                for ( size_t i = 0, e_i = this->diag_A_.size ( ); i != e_i; ++i )
                {
                    int diag_num = offset_A + i;
                    this->A_.GetValues ( &diag_num, 1, &diag_num, 1, &( this->diag_A_[i] ) );
                }
            }
            //*****************************************************************
            // Initialize diag_D_
            //*****************************************************************
            if ( block_id == 3 )
            {
                this->diag_D_.clear ( );
                this->diag_D_.resize ( this->sparsity_D_.size ( ), 0. );

                const int offset_D = this->offsets_block_two_[this->my_rank_];
                for ( size_t i = 0, e_i = diag_D_.size ( ); i != e_i; ++i )
                {
                    int diag_num = offset_D + i;
                    this->D_.GetValues ( &diag_num, 1, &diag_num, 1, &( this->diag_D_[i] ) );
                }
            }

            //*****************************************************************
            // Ensure consistency with general preconditioner data structure 
            //*****************************************************************
            this->op_ = &op;
            this->SetModifiedOperator ( true );
            this->SetUseSolverOperator ( false );
        }

        template<class LAD>
        void TriBlock<LAD>::SetupOperatorE ( OperatorType& op )
        {
            if ( !this->E_is_initialized_ )
            {
                std::vector<int> rows_diag, cols_diag, rows_offdiag, cols_offdiag;
                for ( int i = 0; i < this->sparsity_D_.size ( ); ++i )
                {
                    for ( int j = 0; j < this->sparsity_D_[i].size ( ); ++j )
                    {
                        if (
                             ( this->sparsity_D_[i][j] >= this->offsets_block_two_[this->my_rank_] )
                             && ( this->sparsity_D_[i][j] < this->offsets_block_two_[this->my_rank_ + 1] )
                             )
                        {
                            rows_diag.push_back ( this->offsets_block_two_[this->my_rank_] + i );
                            cols_diag.push_back ( this->sparsity_D_[i][j] );
                        }
                        else
                        {
                            rows_offdiag.push_back ( this->offsets_block_two_[this->my_rank_] + i );
                            cols_offdiag.push_back ( this->sparsity_D_[i][j] );
                        }
                    }
                }

                E_.InitStructure (
                                   vec2ptr ( rows_diag ),
                                   vec2ptr ( cols_diag ),
                                   cols_diag.size ( ),
                                   vec2ptr ( rows_offdiag ),
                                   vec2ptr ( cols_offdiag ),
                                   cols_offdiag.size ( )
                                   );
                this->E_is_initialized_ = true;
            }

            //Q_.Zeros ( );
            int row_index_block = 0;
            for ( typename std::map<int, SortedArray<int> >::const_iterator it = this->couplings_D_.begin ( ),
                  e_it = this->couplings_D_.end ( );
                  it != e_it;
                  ++it )
            {
                // Get values from system matrix op
                int row_index_system = it->first;
                std::vector<ValueType> vals ( it->second.size ( ), 0. );
                op.GetValues ( &row_index_system, 1, vec2ptr ( it->second.data ( ) ), it->second.size ( ), vec2ptr ( vals ) );

                // Set values in block matrix Q
                int block_index_system = row_index_block + this->offsets_block_two_[this->my_rank_];
                this->E_.SetValues ( &block_index_system, 1, vec2ptr ( this->sparsity_D_[row_index_block] ), it->second.size ( ), vec2ptr ( vals ) );
                ++row_index_block;
            }
            this->E_is_set_ = true;
            this->E_modified_ = true;
            this->E_passed2solver_ = false;

            //*****************************************************************
            // Ensure consistency with general preconditioner data structure 
            //*****************************************************************
            this->SetState ( false );
        }

        template<class LAD>
        void TriBlock<LAD>::SetupSolverE ( LinearSolver<LAD>& solver )
        {
            this->solver_E_ = &solver;
        }

        template<class LAD>
        void TriBlock<LAD>::FilterSubVector ( VectorType* y, VectorType* rhs_temp, std::vector<int>* indexset_two )
        {
            std::vector<ValueType> val_temp;
            rhs_temp->Zeros ( );
            val_temp.resize ( this->mapb2s_two_.size ( ), 0. );
            y->GetValues ( vec2ptr ( *indexset_two ), indexset_two->size ( ), vec2ptr ( val_temp ) );
            rhs_temp->SetValues ( vec2ptr ( this->mapb2s_two_ ), this->mapb2s_two_.size ( ), vec2ptr ( val_temp ) );

            rhs_temp->Update ( );
            this->non_lin_op_->ApplyFilter ( *rhs_temp );

            val_temp.clear ( );
            val_temp.resize ( this->mapb2s_two_.size ( ), 0. );
            rhs_temp->GetValues ( vec2ptr ( this->mapb2s_two_ ), this->mapb2s_two_.size ( ), vec2ptr ( val_temp ) );
            y->SetValues ( vec2ptr ( *indexset_two ), indexset_two->size ( ), vec2ptr ( val_temp ) );
        }

        template<class LAD>
        void TriBlock<LAD>::FilterVector ( VectorType* y )
        {
            this->non_lin_op_->ApplyFilter ( *y );
        }

        template<class LAD>
        void TriBlock<LAD>::Build ( )
        {
            if ( this->print_level_ > 2 )
            {
                LOG_INFO ( "[" << this->nested_level_ << "] Build Solver", 1 );
            }

            if ( this->solver_A_ != NULL )
            {
                if ( this->block_op_[0] == 1 )
                {
                    if ( this->A_modified_ || !( this->A_passed2solver_ ) )
                    {
                        if ( this->solver_A_->GetPreconditioner ( ) != NULL )
                        {
                            this->solver_A_->GetPreconditioner ( )->SetupOperator ( this->A_ );
                        }
                        this->solver_A_->SetupOperator ( this->A_ );
                        this->solver_A_->Build ( );

                        this->A_modified_ = false;
                        this->A_passed2solver_ = true;
                    }
                }
            }
            else
            {
                LOG_ERROR ( "[" << this->nested_level_ << "] No solver for submatrix A defined!!!" );
                exit ( -1 );
            }
            if ( this->use_press_conv_diff_ )
            {
                if ( this->scaling_S_ != 0.0 )
                {
                    if ( this->solver_Q_ != NULL )
                    {
                        if ( this->Q_modified_ || !( this->Q_passed2solver_ ) )
                        {
                            if ( this->solver_Q_->GetPreconditioner ( ) != NULL )
                            {
                                this->solver_Q_->GetPreconditioner ( )->SetupOperator ( this->Q_ );
                            }
                            this->solver_Q_->SetupOperator ( this->Q_ );
                            this->solver_Q_->Build ( );
                            this->Q_modified_ = false;
                            this->Q_passed2solver_ = true;
                        }
                    }
                    else
                    {
                        LOG_ERROR ( "[" << this->nested_level_ << "] No solver for operator Q defined!!!" );
                        exit ( -1 );
                    }

                    if ( this->solver_H_ != NULL )
                    {
                        if ( this->H_modified_ || !( this->H_passed2solver_ ) )
                        {
                            if ( this->solver_H_->GetPreconditioner ( ) != NULL )
                            {
                                this->solver_H_->GetPreconditioner ( )->SetupOperator ( this->H_ );
                            }
                            this->solver_H_->SetupOperator ( this->H_ );
                            this->solver_H_->Build ( );
                            this->H_modified_ = false;
                            this->H_passed2solver_ = true;
                        }
                    }
                    else
                    {
                        LOG_ERROR ( "[" << this->nested_level_ << "] No solver for operator H defined !!!" );
                        exit ( -1 );
                    }

                    if ( !this->Q_is_set_ )
                    {
                        LOG_ERROR ( "[" << this->nested_level_ << "] No operator Q defined!!!" );
                        exit ( -1 );
                    }
                    if ( !this->F_is_set_ )
                    {
                        LOG_ERROR ( "[" << this->nested_level_ << "] No operator F defined!!!" );
                        exit ( -1 );
                    }
                    if ( !this->H_is_set_ )
                    {
                        LOG_ERROR ( "[" << this->nested_level_ << "] No operator H defined!!!" );
                        exit ( -1 );
                    }
                }
                if ( this->scaling_E_ != 0.0 )
                {
                    if ( this->solver_E_ != NULL )
                    {
                        if ( this->E_modified_ || !( this->E_passed2solver_ ) )
                        {
                            if ( this->solver_E_->GetPreconditioner ( ) != NULL )
                            {
                                this->solver_E_->GetPreconditioner ( )->SetupOperator ( this->D_ );
                            }
                            this->solver_E_->SetupOperator ( this->E_ );
                            this->solver_E_->Build ( );
                            this->E_modified_ = false;
                            this->E_passed2solver_ = true;
                        }
                    }
                    else
                    {
                        LOG_ERROR ( "[" << this->nested_level_ << "] No solver for operator E defined !!!" );
                        exit ( -1 );
                    }
                }
            }
            else
            {
                if ( this->solver_D_ != NULL )
                {
                    if ( this->block_op_[3] == 1 )
                    {
                        if ( this->D_modified_ || !( this->D_passed2solver_ ) )
                        {
                            if ( this->solver_D_->GetPreconditioner ( ) != NULL )
                            {
                                this->solver_D_->GetPreconditioner ( )->SetupOperator ( this->D_ );
                            }
                            this->solver_D_->SetupOperator ( this->D_ );
                            this->solver_D_->Build ( );
                            this->D_modified_ = false;
                            this->D_passed2solver_ = true;
                        }
                    }
                }
            }
            this->SetModifiedOperator ( false );
            this->SetState ( true );
        }

        template<class LAD>
        LinearSolverState TriBlock<LAD>::Solve ( const VectorType& b, VectorType* x )
        {
            if ( !this->GetState ( ) )
            {
                this->Build ( );
            }

            this->x_.Zeros ( );
            this->y_.Zeros ( );

            //*****************************************************************
            // Create block index sets
            //*****************************************************************
            std::vector<int> indexset_one ( this->mapb2s_one_.size ( ) );
            for ( int i = 0, e_i = indexset_one.size ( ); i != e_i; ++i )
            {
                indexset_one[i] = i + this->offsets_block_one_[this->my_rank_];
            }

            std::vector<int>indexset_two ( this->mapb2s_two_.size ( ) );
            for ( int i = 0, e_i = indexset_two.size ( ); i != e_i; ++i )
            {
                indexset_two[i] = i + this->offsets_block_two_[this->my_rank_];
            }

            //*****************************************************************
            // Get block parts of vector b
            //*****************************************************************
            std::vector<ValueType> val_temp;
            VectorType rhs_temp, sol_temp;
            rhs_temp.CloneFromWithoutContent ( b );
            sol_temp.CloneFromWithoutContent ( b );

            //*****************************************************************
            // Auxiliary vectors
            //*****************************************************************
            // In block one space
            VectorType r_x;
            r_x.CloneFromWithoutContent ( this->x_ );

            // In block two space
            VectorType r_y;
            r_y.CloneFromWithoutContent ( this->y_ );

            //*****************************************************************
            // Solve single block
            Timer timer;

            // ****************************************************************	
            // Solve Lower Triangular system 
            // | A   0 | * |x|  = |f|
            // | C   D |   |y|    |g|
            // ****************************************************************	
            //int rank;
            //MPI_Comm_rank(MPI_COMM_WORLD, &rank);

            if ( block_op_[1] == 0 && block_op_[2] > 0 )
            {
                if ( this->print_level_ > 1 ) LOG_INFO ( "[" << this->nested_level_ << "] Solve lower block triangular system                        ", 1 );

                // ****************************************************************	
                // solve   x = A^{-1} f
                // ****************************************************************

                if ( block_op_[0] == 0 )
                {
                    // operator of solver A is full matrix -> give global rhs to solver A with zeros in block two 	
                    rhs_temp.Zeros ( );
                    sol_temp.Zeros ( );
                    val_temp.clear ( );

                    val_temp.resize ( this->mapb2s_one_.size ( ), 0. );
                    b.GetValues ( vec2ptr ( this->mapb2s_one_ ), this->mapb2s_one_.size ( ), vec2ptr ( val_temp ) );
                    rhs_temp.SetValues ( vec2ptr ( this->mapb2s_one_ ), this->mapb2s_one_.size ( ), vec2ptr ( val_temp ) );

                    this->SolveBlock ( rhs_temp, sol_temp, 0 );

                    val_temp.clear ( );
                    val_temp.resize ( this->mapb2s_one_.size ( ), 0. );
                    sol_temp.GetValues ( vec2ptr ( this->mapb2s_one_ ), this->mapb2s_one_.size ( ), vec2ptr ( val_temp ) );
                    x->SetValues ( vec2ptr ( this->mapb2s_one_ ), this->mapb2s_one_.size ( ), vec2ptr ( val_temp ) );
                    this->x_.SetValues ( vec2ptr ( indexset_one ), indexset_one.size ( ), vec2ptr ( val_temp ) );

                }
                else
                {
                    // operator of A is submatrix A	-> extract block one values and store in vector f_ 
                    val_temp.resize ( this->mapb2s_one_.size ( ), 0. );
                    b.GetValues ( vec2ptr ( this->mapb2s_one_ ), this->mapb2s_one_.size ( ), vec2ptr ( val_temp ) );
                    this->f_.SetValues ( vec2ptr ( indexset_one ), indexset_one.size ( ), vec2ptr ( val_temp ) );

                    this->SolveBlock ( this->f_, this->x_, 0 );

                    val_temp.clear ( );
                    val_temp.resize ( this->mapb2s_one_.size ( ), 0. );
                    this->x_.GetValues ( vec2ptr ( indexset_one ), indexset_one.size ( ), vec2ptr ( val_temp ) );
                    x->SetValues ( vec2ptr ( this->mapb2s_one_ ), this->mapb2s_one_.size ( ), vec2ptr ( val_temp ) );
                }

                // ****************************************************************	
                // compute g <- g - Cx
                // ****************************************************************	

                val_temp.clear ( );
                val_temp.resize ( this->mapb2s_two_.size ( ), 0. );
                b.GetValues ( vec2ptr ( this->mapb2s_two_ ), this->mapb2s_two_.size ( ), vec2ptr ( val_temp ) );
                this->g_.SetValues ( vec2ptr ( indexset_two ), indexset_two.size ( ), vec2ptr ( val_temp ) );

                timer.start ( );
                this->C_.VectorMult ( this->x_, &r_y );
                this->g_.Axpy ( r_y, static_cast < ValueType > ( -1. ) );
                timer.stop ( );

                if ( this->print_level_ > 1 ) LOG_INFO ( "[" << this->nested_level_ << "] residual update: CPU time                                  ", timer.get_duration ( ) );

                // ****************************************************************	
                // solve   y = D^{-1} g
                // ****************************************************************

                if ( block_op_[3] == 0 && this->use_press_conv_diff_ == false )
                {

                    // operator of solver D is full matrix -> give global rhs to solver D with zeros in block one 	 	
                    rhs_temp.Zeros ( );
                    sol_temp.Zeros ( );
                    val_temp.clear ( );

                    val_temp.resize ( this->mapb2s_two_.size ( ), 0. );
                    this->g_.GetValues ( vec2ptr ( indexset_two ), indexset_two.size ( ), vec2ptr ( val_temp ) );
                    rhs_temp.SetValues ( vec2ptr ( this->mapb2s_two_ ), this->mapb2s_two_.size ( ), vec2ptr ( val_temp ) );

                    this->SolveBlock ( rhs_temp, sol_temp, 3 );

                    val_temp.clear ( );
                    val_temp.resize ( this->mapb2s_two_.size ( ), 0. );
                    sol_temp.GetValues ( vec2ptr ( this->mapb2s_two_ ), this->mapb2s_two_.size ( ), vec2ptr ( val_temp ) );
                    x->SetValues ( vec2ptr ( this->mapb2s_two_ ), this->mapb2s_two_.size ( ), vec2ptr ( val_temp ) );
                }
                else
                {

                    // operator of D is submatrix D	
                    this->SolveBlock ( this->g_, this->y_, 3 );

                    val_temp.clear ( );
                    val_temp.resize ( this->mapb2s_two_.size ( ), 0. );
                    this->y_.GetValues ( vec2ptr ( indexset_two ), indexset_two.size ( ), vec2ptr ( val_temp ) );
                    x->SetValues ( vec2ptr ( this->mapb2s_two_ ), this->mapb2s_two_.size ( ), vec2ptr ( val_temp ) );
                }
            }

            // ****************************************************************	
            // Solve Upper Triangular system 
            // | A   B | * |x|  = |f|
            // | 0   D |   |y|  = |g|
            // ****************************************************************	

            if ( block_op_[1] > 0 && block_op_[2] == 0 )
            {
                if ( this->print_level_ > 1 ) LOG_INFO ( "[" << this->nested_level_ << "] Solve upper block triangular system                        ", 1 );

                // ****************************************************************	
                // solve   y = D^{-1} g
                // ****************************************************************

                if ( block_op_[3] == 0 && this->use_press_conv_diff_ == false )
                {
                    // operator of solver D is full matrix -> give global rhs to solver D with zeros in block one 	
                    rhs_temp.Zeros ( );
                    sol_temp.Zeros ( );
                    val_temp.clear ( );

                    val_temp.resize ( this->mapb2s_two_.size ( ), 0. );
                    b.GetValues ( vec2ptr ( this->mapb2s_two_ ), this->mapb2s_two_.size ( ), vec2ptr ( val_temp ) );
                    rhs_temp.SetValues ( vec2ptr ( this->mapb2s_two_ ), this->mapb2s_two_.size ( ), vec2ptr ( val_temp ) );

                    this->SolveBlock ( rhs_temp, sol_temp, 3 );

                    val_temp.clear ( );
                    val_temp.resize ( this->mapb2s_two_.size ( ), 0. );
                    sol_temp.GetValues ( vec2ptr ( this->mapb2s_two_ ), this->mapb2s_two_.size ( ), vec2ptr ( val_temp ) );
                    x->SetValues ( vec2ptr ( this->mapb2s_two_ ), this->mapb2s_two_.size ( ), vec2ptr ( val_temp ) );
                    this->y_.SetValues ( vec2ptr ( indexset_two ), indexset_two.size ( ), vec2ptr ( val_temp ) );

                }
                else
                {
                    // operator of D is submatrix D	-> extract block two values and store in vector g_ 
                    val_temp.resize ( this->mapb2s_two_.size ( ), 0. );
                    b.GetValues ( vec2ptr ( this->mapb2s_two_ ), this->mapb2s_two_.size ( ), vec2ptr ( val_temp ) );
                    this->g_.SetValues ( vec2ptr ( indexset_two ), indexset_two.size ( ), vec2ptr ( val_temp ) );

                    this->SolveBlock ( this->g_, this->y_, 3 );

                    val_temp.clear ( );
                    val_temp.resize ( this->mapb2s_two_.size ( ), 0. );
                    this->y_.GetValues ( vec2ptr ( indexset_two ), indexset_two.size ( ), vec2ptr ( val_temp ) );
                    x->SetValues ( vec2ptr ( this->mapb2s_two_ ), this->mapb2s_two_.size ( ), vec2ptr ( val_temp ) );
                }

                // ****************************************************************	
                // compute f <- f - By
                // ****************************************************************	

                val_temp.clear ( );
                val_temp.resize ( this->mapb2s_one_.size ( ), 0. );
                b.GetValues ( vec2ptr ( this->mapb2s_one_ ), this->mapb2s_one_.size ( ), vec2ptr ( val_temp ) );
                this->f_.SetValues ( vec2ptr ( indexset_one ), indexset_one.size ( ), vec2ptr ( val_temp ) );

                timer.start ( );
                this->B_.VectorMult ( this->y_, &r_x );
                this->f_.Axpy ( r_x, static_cast < ValueType > ( -1. ) );
                timer.stop ( );

                if ( this->print_level_ > 1 ) LOG_INFO ( "[" << this->nested_level_ << "] residual update: CPU time                                  ", timer.get_duration ( ) );

                // ****************************************************************	
                // solve   x = A^{-1} f
                // ****************************************************************

                if ( block_op_[0] == 0 )
                {
                    // operator of solver A is full matrix -> give global rhs to solver A with zeros in block two 	 	
                    rhs_temp.Zeros ( );
                    sol_temp.Zeros ( );
                    val_temp.clear ( );

                    val_temp.resize ( this->mapb2s_one_.size ( ), 0. );
                    this->f_.GetValues ( vec2ptr ( indexset_one ), indexset_one.size ( ), vec2ptr ( val_temp ) );
                    rhs_temp.SetValues ( vec2ptr ( this->mapb2s_one_ ), this->mapb2s_one_.size ( ), vec2ptr ( val_temp ) );

                    this->SolveBlock ( rhs_temp, sol_temp, 0 );

                    val_temp.clear ( );
                    val_temp.resize ( this->mapb2s_one_.size ( ), 0. );
                    sol_temp.GetValues ( vec2ptr ( this->mapb2s_one_ ), this->mapb2s_one_.size ( ), vec2ptr ( val_temp ) );
                    x->SetValues ( vec2ptr ( this->mapb2s_one_ ), this->mapb2s_one_.size ( ), vec2ptr ( val_temp ) );
                }
                else
                {
                    // operator of A is submatrix A	
                    this->SolveBlock ( this->f_, this->x_, 0 );

                    val_temp.clear ( );
                    val_temp.resize ( this->mapb2s_one_.size ( ), 0. );
                    this->x_.GetValues ( vec2ptr ( indexset_one ), indexset_one.size ( ), vec2ptr ( val_temp ) );
                    x->SetValues ( vec2ptr ( this->mapb2s_one_ ), this->mapb2s_one_.size ( ), vec2ptr ( val_temp ) );
                }
            }

            // ****************************************************************	
            // Solve Diagonal system 
            // | A   0 | * |x|  = |f|
            // | 0   D |   |y|  = |g|
            // ****************************************************************	

            if ( this->block_op_[1] == 0 && this->block_op_[2] == 0 )
            {
                if ( this->print_level_ > 1 ) LOG_INFO ( "[" << this->nested_level_ << "] Solve block diagonal system                                ", 1 );

                // ****************************************************************	
                // solve   x = A^{-1} f
                // ****************************************************************

                if ( block_op_[0] == 0 )
                {
                    // operator of solver A is full matrix -> give global rhs to solver A with zeros in block two 	
                    rhs_temp.Zeros ( );
                    sol_temp.Zeros ( );
                    val_temp.clear ( );

                    val_temp.resize ( this->mapb2s_one_.size ( ), 0. );
                    b.GetValues ( vec2ptr ( this->mapb2s_one_ ), this->mapb2s_one_.size ( ), vec2ptr ( val_temp ) );
                    rhs_temp.SetValues ( vec2ptr ( this->mapb2s_one_ ), this->mapb2s_one_.size ( ), vec2ptr ( val_temp ) );

                    this->SolveBlock ( rhs_temp, sol_temp, 0 );

                    val_temp.clear ( );
                    val_temp.resize ( this->mapb2s_one_.size ( ), 0. );
                    sol_temp.GetValues ( vec2ptr ( this->mapb2s_one_ ), this->mapb2s_one_.size ( ), vec2ptr ( val_temp ) );
                    x->SetValues ( vec2ptr ( this->mapb2s_one_ ), this->mapb2s_one_.size ( ), vec2ptr ( val_temp ) );
                }
                else
                {
                    // operator of A is submatrix A	-> extract block one values and store in vector f_ 
                    val_temp.resize ( this->mapb2s_one_.size ( ), 0. );
                    b.GetValues ( vec2ptr ( this->mapb2s_one_ ), this->mapb2s_one_.size ( ), vec2ptr ( val_temp ) );
                    this->f_.SetValues ( vec2ptr ( indexset_one ), indexset_one.size ( ), vec2ptr ( val_temp ) );

                    this->SolveBlock ( this->f_, this->x_, 0 );

                    val_temp.clear ( );
                    val_temp.resize ( this->mapb2s_one_.size ( ), 0. );
                    this->x_.GetValues ( vec2ptr ( indexset_one ), indexset_one.size ( ), vec2ptr ( val_temp ) );
                    x->SetValues ( vec2ptr ( this->mapb2s_one_ ), this->mapb2s_one_.size ( ), vec2ptr ( val_temp ) );
                }

                // ****************************************************************	
                // solve   y = D^{-1} g
                // ****************************************************************

                if ( block_op_[3] == 0 && this->use_press_conv_diff_ == false )
                {
                    // operator of solver D is full matrix -> give global rhs to solver D with zeros in block one 	 	
                    rhs_temp.Zeros ( );
                    sol_temp.Zeros ( );
                    val_temp.clear ( );

                    val_temp.resize ( this->mapb2s_two_.size ( ), 0. );
                    b.GetValues ( vec2ptr ( this->mapb2s_two_ ), this->mapb2s_two_.size ( ), vec2ptr ( val_temp ) );
                    rhs_temp.SetValues ( vec2ptr ( this->mapb2s_two_ ), this->mapb2s_two_.size ( ), vec2ptr ( val_temp ) );

                    this->SolveBlock ( rhs_temp, sol_temp, 3 );

                    val_temp.clear ( );
                    val_temp.resize ( this->mapb2s_two_.size ( ), 0. );
                    sol_temp.GetValues ( vec2ptr ( this->mapb2s_two_ ), this->mapb2s_two_.size ( ), vec2ptr ( val_temp ) );
                    x->SetValues ( vec2ptr ( this->mapb2s_two_ ), this->mapb2s_two_.size ( ), vec2ptr ( val_temp ) );
                }
                else
                {
                    // operator of D is submatrix D	
                    val_temp.resize ( this->mapb2s_two_.size ( ), 0. );
                    b.GetValues ( vec2ptr ( this->mapb2s_two_ ), this->mapb2s_two_.size ( ), vec2ptr ( val_temp ) );
                    this->g_.SetValues ( vec2ptr ( indexset_two ), indexset_two.size ( ), vec2ptr ( val_temp ) );

                    this->SolveBlock ( this->g_, this->y_, 3 );

                    val_temp.clear ( );
                    val_temp.resize ( this->mapb2s_two_.size ( ), 0. );
                    this->y_.GetValues ( vec2ptr ( indexset_two ), indexset_two.size ( ), vec2ptr ( val_temp ) );
                    x->SetValues ( vec2ptr ( this->mapb2s_two_ ), this->mapb2s_two_.size ( ), vec2ptr ( val_temp ) );
                }
            }

            // ****************************************************************	
            // Can't Solve Full system 
            // | A   B | * |x|  = |f|
            // | C   D |   |y|  = |g|
            // ****************************************************************	

            if ( block_op_[1] > 0 && block_op_[2] > 0 )
            {
                LOG_ERROR ( "[" << this->nested_level_ << "] Matrix is neither block triangular nor block diagonal, idiot!!!" );
                exit ( -1 );
            }

            // Filter solution
            if ( this->filter_solution_ ) FilterVector ( x );

            // Cleanup
            indexset_one.clear ( );
            indexset_two.clear ( );
            val_temp.clear ( );

            r_x.Clear ( );
            r_y.Clear ( );
        }

        template<class LAD>
        void TriBlock<LAD>::SolveBlock ( const VectorType &b, VectorType &x, const int block_id )
        {
            if ( this->is_simple_precond_ )
            {
                if ( block_id == 0 )
                {
                    std::vector<ValueType> b_loc ( this->diag_A_.size ( ), 0. );
                    std::vector<ValueType> x_loc ( this->diag_A_.size ( ), 0. );

                    const int offset_A = this->offsets_block_one_[this->my_rank_];

                    std::vector<int> indices_A ( this->diag_A_.size ( ), 0 );
                    for ( size_t i = 0, e_i = this->diag_A_.size ( ); i != e_i; ++i )
                    {
                        indices_A[i] = offset_A + i;
                    }

                    b.GetValues ( vec2ptr ( indices_A ), indices_A.size ( ), vec2ptr ( b_loc ) );

                    for ( size_t i = 0, e_i = this->diag_A_.size ( ); i != e_i; ++i )
                    {
                        x_loc[i] = b_loc[i] / this->diag_A_[i];
                    }

                    x.SetValues ( vec2ptr ( indices_A ), indices_A.size ( ), vec2ptr ( x_loc ) );

                    b_loc.clear ( );
                    x_loc.clear ( );
                }
                if ( block_id == 3 )
                {
                    std::vector<ValueType> b_loc ( this->diag_D_.size ( ), 0. );
                    std::vector<ValueType> x_loc ( this->diag_D_.size ( ), 0. );

                    const int offset_D = this->offsets_block_one_[this->my_rank_];

                    std::vector<int> indices_D ( this->diag_D_.size ( ), 0 );
                    for ( size_t i = 0, e_i = this->diag_D_.size ( ); i != e_i; ++i )
                    {
                        indices_D[i] = offset_D + i;
                    }

                    b.GetValues ( vec2ptr ( indices_D ), indices_D.size ( ), vec2ptr ( b_loc ) );

                    for ( size_t i = 0, e_i = this->diag_D_.size ( ); i != e_i; ++i )
                    {
                        x_loc[i] = b_loc[i] / this->diag_D_[i];
                    }

                    x.SetValues ( vec2ptr ( indices_D ), indices_D.size ( ), vec2ptr ( x_loc ) );
                    b_loc.clear ( );
                    x_loc.clear ( );
                }
            }
            else
            {
                Timer timer;
                if ( block_id == 0 )
                {
                    timer.start ( );
                    this->solver_A_->Solve ( b, &x );
                    timer.stop ( );
                    if ( this->print_level_ == 1 )
                    {
                        LOG_INFO ( "[" << this->nested_level_ << "] block A solver : CPU time                                  ", timer.get_duration ( ) );
                    }
                    if ( this->print_level_ == 2 )
                    {
                        LOG_INFO ( "[" << this->nested_level_ << "] block A solver : CPU time                                  ", timer.get_duration ( ) );
                        LOG_INFO ( "[" << this->nested_level_ << "] block A solver : Iter                                      ", this->solver_A_->iter ( ) );
                        LOG_INFO ( "[" << this->nested_level_ << "] block A solver : final rel. res.                           ", this->solver_A_->res ( ) );
                    }
                }
                if ( block_id == 3 )
                {
                    timer.start ( );
                    if ( this->use_press_conv_diff_ )
                    {
                        this->ApplyPressureConvDiff ( b, &x );
                        timer.stop ( );
                        if ( this->print_level_ == 1 )
                        {
                            LOG_INFO ( "[" << this->nested_level_ << "] approx Schur   : CPU time                                  ", timer.get_duration ( ) );
                        }
                        if ( this->print_level_ == 2 )
                        {
                            LOG_INFO ( "[" << this->nested_level_ << "] approx Schur   : CPU time                                  ", timer.get_duration ( ) );
                            if ( this->solver_E_ != NULL ) LOG_INFO ( "[" << this->nested_level_ << "] block E solver : Iter                                      ", this->solver_E_->iter ( ) );
                            if ( this->solver_E_ != NULL ) LOG_INFO ( "[" << this->nested_level_ << "] block E solver : final rel. res.                           ", this->solver_E_->res ( ) );
                            if ( this->solver_Q_ != NULL ) LOG_INFO ( "[" << this->nested_level_ << "] block Q solver : Iter                                      ", this->solver_Q_->iter ( ) );
                            if ( this->solver_Q_ != NULL ) LOG_INFO ( "[" << this->nested_level_ << "] block Q solver : final rel. res.                           ", this->solver_Q_->res ( ) );
                            if ( this->solver_H_ != NULL ) LOG_INFO ( "[" << this->nested_level_ << "] block H solver : Iter                                      ", this->solver_H_->iter ( ) );
                            if ( this->solver_H_ != NULL ) LOG_INFO ( "[" << this->nested_level_ << "] block H solver : final rel. res.                           ", this->solver_H_->res ( ) );
                        }
                    }
                    else
                    {
                        this->solver_D_->Solve ( b, &x );
                        timer.stop ( );
                        if ( this->print_level_ == 1 )
                        {
                            LOG_INFO ( "[" << this->nested_level_ << "] block D solver : CPU time                                  ", timer.get_duration ( ) );
                        }
                        if ( this->print_level_ == 2 )
                        {
                            LOG_INFO ( "[" << this->nested_level_ << "] block D solver : CPU time                                  ", timer.get_duration ( ) );
                            LOG_INFO ( "[" << this->nested_level_ << "] block D solver : Iter                                      ", this->solver_D_->iter ( ) );
                            LOG_INFO ( "[" << this->nested_level_ << "] block D solver : final rel. res.                           ", this->solver_D_->res ( ) );
                        }
                    }
                }
            }
        }

        template<class LAD>
        void TriBlock<LAD>::ApplyPressureConvDiff ( const VectorType &b, VectorType *x )
        {
            VectorType y, z;
            x->Zeros ( );
            y.CloneFromWithoutContent ( *x );
            z.CloneFromWithoutContent ( *x );

            // compute -scaling_S_ * Q^{-1} * F * H^{-1}b
            if ( this->scaling_S_ != 0.0 )
            {
                // y = 	H^{-1}b
                this->solver_H_->Solve ( b, &y );

                // z = F * y	
                this->F_.VectorMult ( y, &z );

                // x = Q^{-1} z	
                this->solver_Q_->Solve ( z, x );

                // x = -scaling_S * x
                x->Scale ( static_cast < ValueType > ( -1. ) * this->scaling_S_ );
            }

            // compute x = scaling_E * E^{-1}b + x	
            if ( this->scaling_E_ != 0.0 && this->solver_E_ != NULL )
            {
                // y = E^{-1} b 	
                y.Zeros ( );
                this->solver_E_->Solve ( b, &y );

                // x = x + scaling_E * y
                x->Axpy ( y, static_cast < ValueType > ( 1. ) * this->scaling_E_ );
            }
            z.Clear ( );
            y.Clear ( );
        }

        // template instantiation
        template class TriBlock<LADescriptorHypreD>;
    }
}
