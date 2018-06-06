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

#include "dof/dof_partition.h"
#include "mesh/iterator.h"
#include "common/macros.h"
#include "common/log.h"
#include <cassert>
#include <vector>
#include <map>

const int DEBUG_LEVEL = 0;

/// \author Michael Schick<br>Martin Baumann<br>Simon Gawlok

namespace hiflow
{
    namespace doffem
    {

        template<class DataType>
        DofPartition<DataType>::DofPartition ( ) : DegreeOfFreedom<DataType>::DegreeOfFreedom ( )
        {
            my_dof_offset_ = 0;

            // Default Value is COMM WORLD
            comm_ = MPI_COMM_WORLD;
            MPI_Comm_size ( comm_, &number_of_subdomains_ );
            ndofs_on_sd_.resize ( number_of_subdomains_, -1 );
            MPI_Comm_rank ( comm_, &my_subdomain_ );
            shared_subdomains_.resize ( number_of_subdomains_, false );
        }

        template<class DataType>
        DofPartition<DataType>::~DofPartition ( )
        {
        }

        template<class DataType>
        void DofPartition<DataType>::set_mpi_comm ( const MPI_Comm& comm )
        {
            comm_ = comm;
            MPI_Comm_size ( comm, &number_of_subdomains_ );

            ndofs_on_sd_.resize ( number_of_subdomains_, -1 );
            MPI_Comm_rank ( comm, &my_subdomain_ );
            shared_subdomains_.resize ( number_of_subdomains_, false );
        }

        template<class DataType>
        void DofPartition<DataType>::number ( DOF_ORDERING order )
        {
            DegreeOfFreedom<DataType>::number ( order );

            // Check if sequential case or truely parallel case is used
            if ( number_of_subdomains_ == 1 )
            {
                ndofs_on_sd_[0] = this->get_nb_dofs ( );
            }
            else
            {
                create_ownerships ( );
                renumber ( );

                this->number_of_dofs_total_ = ndofs_on_sd_[my_subdomain_];
            }

            // Calculate number of dofs on the global domain
            ndofs_global_ = 0;
            for ( size_t s = 0, e_s = number_of_subdomains_; s != e_s; ++s )
                ndofs_global_ += ndofs_on_sd_[s];

            // Some essential outputs
            if ( my_subdomain_ == 0 )
            {
                LOG_INFO ( "#Dofs on global domain", ndofs_global_ );
                for ( size_t s = 0, e_s = number_of_subdomains_; s != e_s; ++s )
                {
                    LOG_INFO ( "#Dofs on subdomain " << s, ndofs_on_sd_[s] );
                }
            }
        }

        template<class DataType>
        void DofPartition<DataType>::create_ownerships ( )
        {
            ownership_.clear ( );
            ownership_.resize ( this->get_nb_dofs ( ), -1 );

            for ( mesh::EntityIterator it = this->mesh_->begin ( this->tdim_ ),
                  e_it = this->mesh_->end ( this->tdim_ );
                  it != e_it;
                  ++it )
            {
                int subdomain_of_entity;
                it->get ( "_sub_domain_", &subdomain_of_entity );

                if ( subdomain_of_entity != my_subdomain_ )
                {
                    LOG_DEBUG ( 3, "[" << my_subdomain_ << "] remote cell with id " << it->id ( ) );
                    for ( size_t var = 0, e_var = this->nvar_; var != e_var; ++var )
                    {
                        std::vector<DofID> local_dofs;
                        this->get_dofs_on_cell ( var, it->index ( ), local_dofs );

                        for ( size_t i = 0, e_i = local_dofs.size ( ); i != e_i; ++i )
                        {
                            const DofID ld_i = local_dofs[i];
                            if ( subdomain_of_entity < ownership_[ld_i] || ownership_[ld_i] == -1 )
                                ownership_[ld_i] = subdomain_of_entity;
                        }
                    }
                }
                else
                {
                    LOG_DEBUG ( 3, "[" << my_subdomain_ << "] local cell with id " << it->id ( ) );
                    for ( size_t var = 0, e_var = this->nvar_; var != e_var; ++var )
                    {
                        std::vector<DofID> local_dofs;
                        this->get_dofs_on_cell ( var, it->index ( ), local_dofs );

                        for ( size_t i = 0, e_i = local_dofs.size ( ); i != e_i; ++i )
                        {
                            const DofID ld_i = local_dofs[i];
                            if ( my_subdomain_ < ownership_[ld_i] || ownership_[ld_i] == -1 )
                                ownership_[ld_i] = my_subdomain_;
                        }
                    }
                }
            }
        }

        template<class DataType>
        void DofPartition<DataType>::renumber ( )
        {
            // Communicate number of dofs including ghost layer
            std::vector<int> ndofs_with_ghost ( number_of_subdomains_ );
            ndofs_with_ghost[my_subdomain_] = this->get_nb_dofs ( );

            int ndofs_with_ghost_sent = this->get_nb_dofs ( );

            // Store information about number of dofs including ghost layer
            ndofs_incl_ghost_ = this->get_nb_dofs ( );

            MPI_Allgather ( &ndofs_with_ghost_sent,
                            1, MPI_INT,
                            vec2ptr ( ndofs_with_ghost ),
                            1, MPI_INT, comm_ );

            LOG_DEBUG ( 2, "[" << my_subdomain_ << "]: Number of DOFs on each subdomain: " << string_from_range ( ndofs_with_ghost.begin ( ), ndofs_with_ghost.end ( ) ) );
            LOG_DEBUG ( 2, "[" << my_subdomain_ << "]: Ownership: " << string_from_range ( ownership_.begin ( ), ownership_.end ( ) ) );

            // Calculate temporary dof offset
            int tmp_dof_offset = 0;
            for ( size_t s = 0, e_s = my_subdomain_; s != e_s; ++s )
                tmp_dof_offset += ndofs_with_ghost[s];

            // Fill first permutation to create a local consecutive numbering w.r.t. tmp_dof_offset
            std::vector<int> permutation ( ownership_.size ( ) );
            int dof_number = 0;

            for ( size_t i = 0, e_i = ownership_.size ( ); i != e_i; ++i )
            {
                if ( ownership_[i] == my_subdomain_ )
                {
                    permutation[i] = dof_number + tmp_dof_offset;
                    ++dof_number;
                }
            }

            LOG_DEBUG ( 2, "[" << my_subdomain_ << "]: Permutation size " << permutation.size ( ) << ", content: " << string_from_range ( permutation.begin ( ), permutation.end ( ) ) );

            int ghost_number = 0;
            for ( size_t i = 0, e_i = ownership_.size ( ); i != e_i; ++i )
            {
                if ( ownership_[i] != my_subdomain_ )
                {
                    permutation[i] = dof_number + tmp_dof_offset + ghost_number;
                    ++ghost_number;
                }
            }

            LOG_DEBUG ( 2, "[" << my_subdomain_ << "]: Permutation size " << permutation.size ( ) << ", content: " << string_from_range ( permutation.begin ( ), permutation.end ( ) ) );

            this->apply_permutation ( permutation );

            LOG_DEBUG ( 2, " First permutation done " );

            // Calculate number of dofs which belong to my subdomain
            ndofs_on_sd_[my_subdomain_] = dof_number;

            // Gather number of dofs of all subdomains on all processes
            MPI_Allgather ( &dof_number,
                            1, MPI_INT,
                            vec2ptr ( ndofs_on_sd_ ),
                            1, MPI_INT, comm_ );
            LOG_DEBUG ( 2, "[" << my_subdomain_ << "]: ndofs_on_sd " << string_from_range ( ndofs_on_sd_.begin ( ), ndofs_on_sd_.end ( ) ) );

            // Setting up data structure for communication and management of dofs 
            // which belong to ghost cells on current process.
            // The data structure maps subdomain indices to maps of (ghost) cell 
            // indices. For each (ghost) cell, a map of size nvar_ is hold which 
            // itself contains the (global) dof numbers of this cell
            std::map< int, std::map< int, std::map< int, std::vector< int > > > > numer_ghost;

            // Create subdomain/ghost cell structure of numer_ghost
            for ( mesh::EntityIterator it = this->mesh_->begin ( this->tdim_ ),
                  e_it = this->mesh_->end ( this->tdim_ );
                  it != e_it;
                  ++it )
            {
                int subdomain_index;
                it->get ( "_sub_domain_", &subdomain_index );

                if ( subdomain_index != my_subdomain_ )
                {
                    // Access/create (on first access) map for remote subdomain
                    numer_ghost[subdomain_index];
                    int ghost_cell_index;
                    it->get ( "_remote_index_", &ghost_cell_index );
                    // Access/create (on first access) map entry for ghost cell
                    numer_ghost[subdomain_index][ghost_cell_index];
                }
            }

            // Exchange dof numbers and ids
            for ( int i = 0; i < this->number_of_subdomains_; ++i )
            {
                for ( int j = 0; j < this->number_of_subdomains_; ++j )
                {
                    if ( i != j )
                    {
                        // If this process is "sender"
                        if ( i == this->my_subdomain_ )
                        {
                            MPI_Status stat;

                            // Send number of requested ghost cells
                            int num_ghost_cells = 0;
                            if ( numer_ghost.find ( j ) != numer_ghost.end ( ) )
                            {
                                num_ghost_cells = numer_ghost[j].size ( );
                            }
                            MPI_Send ( &num_ghost_cells, 1, MPI_INT, j, 0, this->comm_ );

                            // Only proceed if num_ghost_cells > 0
                            if ( num_ghost_cells > 0 )
                            {
                                // Send ghost cell indices
                                std::vector<int> ghost_indices;
                                for ( typename std::map< int, std::map< int, std::vector< int > > >::const_iterator it = numer_ghost[j].begin ( ),
                                      e_it = numer_ghost[j].end ( );
                                      it != e_it;
                                      ++it )
                                {
                                    ghost_indices.push_back ( it->first );
                                }

                                MPI_Send ( vec2ptr ( ghost_indices ), num_ghost_cells, MPI_INT, j, 1, this->comm_ );

                                // Receive number of dofs for all variables and all ghost cells
                                std::vector<int> num_dofs_ghost ( num_ghost_cells * this->nvar_, -1 );

                                MPI_Recv ( vec2ptr ( num_dofs_ghost ), num_ghost_cells * this->nvar_, MPI_INT, j, 2, this->comm_, &stat );

                                // Prepare final numer_ghost structure
                                for ( int k = 0; k < num_ghost_cells; ++k )
                                {
                                    for ( int l = 0; l < this->nvar_; ++l )
                                    {
                                        numer_ghost[j][ghost_indices[k]][l].resize ( num_dofs_ghost[k * this->nvar_ + l] );
                                    }
                                }

                                //prepare data to receive
                                int num_ind = 0;
                                for ( int k = 0; k < num_ghost_cells; ++k )
                                {
                                    for ( int l = 0; l < this->nvar_; ++l )
                                    {
                                        num_ind += num_dofs_ghost[k * this->nvar_ + l];
                                    }
                                }

                                std::vector<int> recv_indices ( num_ind );

                                MPI_Recv ( vec2ptr ( recv_indices ), num_ind, MPI_INT, j, 3, this->comm_, &stat );

                                // Unpack received data
                                int ind = 0;
                                for ( int k = 0; k < num_ghost_cells; ++k )
                                {
                                    for ( int l = 0; l < this->nvar_; ++l )
                                    {
                                        for ( int m = 0; m < num_dofs_ghost[k * this->nvar_ + l]; ++m )
                                        {
                                            numer_ghost[j][ghost_indices[k]][l][m] = recv_indices[ind];
                                            ++ind;
                                        }
                                    }
                                }
                            }
                        }

                        // If this process is "receiver"
                        if ( j == this->my_subdomain_ )
                        {
                            MPI_Status stat;

                            // Receive number of requested ghost cells
                            int num_ghost_cells = -1;
                            MPI_Recv ( &num_ghost_cells, 1, MPI_INT, i, 0, this->comm_, &stat );

                            // Only proceed if num_ghost_cells > 0
                            if ( num_ghost_cells > 0 )
                            {
                                // Receive ghost cell indices
                                std::vector<int> ghost_indices ( num_ghost_cells, -1 );

                                MPI_Recv ( vec2ptr ( ghost_indices ), num_ghost_cells, MPI_INT, i, 1, this->comm_, &stat );

                                // Send number of dofs for all variables and all ghost cells
                                std::vector<int> num_dofs_ghost ( num_ghost_cells * this->nvar_, -1 );

                                int num_ind = 0;

                                for ( int k = 0; k < num_ghost_cells; ++k )
                                {
                                    for ( int l = 0; l < this->nvar_; ++l )
                                    {
                                        num_dofs_ghost[k * this->nvar_ + l] = this->get_nb_dofs_on_cell ( l, ghost_indices[k] );
                                        num_ind += num_dofs_ghost[k * this->nvar_ + l];
                                    }
                                }

                                MPI_Send ( vec2ptr ( num_dofs_ghost ), num_ghost_cells * this->nvar_, MPI_INT, i, 2, this->comm_ );

                                // Prepare data to send
                                std::vector<int> send_indices;
                                send_indices.reserve ( num_ind );
                                for ( int k = 0; k < num_ghost_cells; ++k )
                                {
                                    for ( int l = 0; l < this->nvar_; ++l )
                                    {
                                        std::vector<int> dof_indices;
                                        this->get_dofs_on_cell ( l, ghost_indices[k], dof_indices );
                                        for ( int m = 0; m < dof_indices.size ( ); ++m )
                                        {
                                            send_indices.push_back ( dof_indices[m] );
                                        }
                                    }
                                }

                                // Send data
                                MPI_Send ( vec2ptr ( send_indices ), num_ind, MPI_INT, i, 3, this->comm_ );
                            }
                        }
                    }
                }
            }

            // First exchange of temporary Dof Ids for ghost layer dofs

            int max_dof_id = *( std::max_element ( this->numer_.begin ( ), this->numer_.end ( ) ) );

            std::vector<int> tmp_permutation ( max_dof_id + 1 );
            for ( size_t i = 0, e_i = tmp_permutation.size ( ); i != e_i; ++i )
                tmp_permutation[i] = i;

            for ( mesh::EntityIterator it = this->mesh_->begin ( this->tdim_ ),
                  e_it = this->mesh_->end ( this->tdim_ );
                  it != e_it;
                  ++it )
            {
                int subdomain_index;
                it->get ( "_sub_domain_", &subdomain_index );

                if ( subdomain_index != my_subdomain_ )
                {
                    int ghost_cell_index;
                    it->get ( "_remote_index_", &ghost_cell_index );

                    int hostile_tmp_dof_offset = 0;
                    for ( size_t s = 0, e_s = subdomain_index; s != e_s; ++s )
                        hostile_tmp_dof_offset += ndofs_with_ghost[s];

                    for ( size_t var = 0, e_var = this->nvar_; var != e_var; ++var )
                    {
                        int size = this->fe_manager_->get_fe_on_cell ( it->index ( ), var )->get_nb_dof_on_cell ( );

                        // Get temporary dof ids from other subdomain
                        std::vector<DofID> ghost_layer_dofs ( size );
                        for ( size_t i = 0, e_i = ghost_layer_dofs.size ( ); i != e_i; ++i )
                        {
                            ghost_layer_dofs[i] = numer_ghost[subdomain_index][ghost_cell_index][var][i];
                        }

                        // Get corresponding dof ids on ghost layer, which need to be updated
                        std::vector<DofID> critical_ghost_layer_dofs;
                        this->get_dofs_on_cell ( var, it->index ( ), critical_ghost_layer_dofs );

                        for ( size_t i = 0, e_i = critical_ghost_layer_dofs.size ( ); i != e_i; ++i )
                        {
                            const int cgld_i = critical_ghost_layer_dofs[i];
                            if ( cgld_i >= tmp_dof_offset + ndofs_on_sd_[my_subdomain_]
                                 || cgld_i < tmp_dof_offset )
                            {
                                const int gld_i = ghost_layer_dofs[i];
                                if ( gld_i >= hostile_tmp_dof_offset
                                     && gld_i < hostile_tmp_dof_offset + ndofs_on_sd_[subdomain_index] )
                                {
                                    assert ( cgld_i >= 0 && cgld_i < tmp_permutation.size ( ) );
                                    tmp_permutation[cgld_i] = gld_i;
                                }
                            }
                        }
                    }
                }
            }

            this->apply_permutation ( tmp_permutation );

            LOG_DEBUG ( 2, " Second permutation done " );

            // Update numer field for all subdomains
            numer_ghost.clear ( );

            // Create subdomain/ghost cell structure of numer_ghost
            for ( mesh::EntityIterator it = this->mesh_->begin ( this->tdim_ ),
                  e_it = this->mesh_->end ( this->tdim_ );
                  it != e_it;
                  ++it )
            {
                int subdomain_index;
                it->get ( "_sub_domain_", &subdomain_index );

                if ( subdomain_index != my_subdomain_ )
                {
                    // Access/create (on first access) map for remote subdomain
                    numer_ghost[subdomain_index];
                    int ghost_cell_index;
                    it->get ( "_remote_index_", &ghost_cell_index );
                    // Access/create (on first access) map entry for ghost cell
                    numer_ghost[subdomain_index][ghost_cell_index];
                }
            }

            // Exchange dof numbers and ids
            for ( int i = 0; i < this->number_of_subdomains_; ++i )
            {
                for ( int j = 0; j < this->number_of_subdomains_; ++j )
                {
                    if ( i != j )
                    {
                        // If this process is "sender"
                        if ( i == this->my_subdomain_ )
                        {
                            MPI_Status stat;

                            // Send number of requested ghost cells
                            int num_ghost_cells = 0;
                            if ( numer_ghost.find ( j ) != numer_ghost.end ( ) )
                            {
                                num_ghost_cells = numer_ghost[j].size ( );
                            }
                            MPI_Send ( &num_ghost_cells, 1, MPI_INT, j, 0, this->comm_ );

                            // Only proceed if num_ghost_cells > 0
                            if ( num_ghost_cells > 0 )
                            {
                                // Send ghost cell indices
                                std::vector<int> ghost_indices;
                                for ( typename std::map< int, std::map< int, std::vector< int > > >::const_iterator it = numer_ghost[j].begin ( ),
                                      e_it = numer_ghost[j].end ( );
                                      it != e_it;
                                      ++it )
                                {
                                    ghost_indices.push_back ( it->first );
                                }

                                MPI_Send ( vec2ptr ( ghost_indices ), num_ghost_cells, MPI_INT, j, 1, this->comm_ );

                                // Receive number of dofs for all variables and all ghost cells
                                std::vector<int> num_dofs_ghost ( num_ghost_cells * this->nvar_, -1 );

                                MPI_Recv ( vec2ptr ( num_dofs_ghost ), num_ghost_cells * this->nvar_, MPI_INT, j, 2, this->comm_, &stat );

                                // Prepare final numer_ghost structure
                                for ( int k = 0; k < num_ghost_cells; ++k )
                                {
                                    for ( int l = 0; l < this->nvar_; ++l )
                                    {
                                        numer_ghost[j][ghost_indices[k]][l].resize ( num_dofs_ghost[k * this->nvar_ + l] );
                                    }
                                }

                                // Receive ghost dof indices
                                for ( int k = 0; k < num_ghost_cells; ++k )
                                {
                                    for ( int l = 0; l < this->nvar_; ++l )
                                    {
                                        MPI_Recv ( vec2ptr ( numer_ghost[j][ghost_indices[k]][l] ), num_dofs_ghost[k * this->nvar_ + l], MPI_INT, j, 0, this->comm_, &stat );
                                    }
                                }
                            }
                        }

                        // If this process is "receiver"
                        if ( j == this->my_subdomain_ )
                        {
                            MPI_Status stat;

                            // Receive number of requested ghost cells
                            int num_ghost_cells = -1;
                            MPI_Recv ( &num_ghost_cells, 1, MPI_INT, i, 0, this->comm_, &stat );

                            // Only proceed if num_ghost_cells > 0
                            if ( num_ghost_cells > 0 )
                            {
                                // Receive ghost cell indices
                                std::vector<int> ghost_indices ( num_ghost_cells, -1 );

                                MPI_Recv ( vec2ptr ( ghost_indices ), num_ghost_cells, MPI_INT, i, 1, this->comm_, &stat );

                                // Send number of dofs for all variables and all ghost cells
                                std::vector<int> num_dofs_ghost ( num_ghost_cells * this->nvar_, -1 );

                                for ( int k = 0; k < num_ghost_cells; ++k )
                                {
                                    for ( int l = 0; l < this->nvar_; ++l )
                                    {
                                        num_dofs_ghost[k * this->nvar_ + l] = this->get_nb_dofs_on_cell ( l, ghost_indices[k] );
                                    }
                                }

                                MPI_Send ( vec2ptr ( num_dofs_ghost ), num_ghost_cells * this->nvar_, MPI_INT, i, 2, this->comm_ );

                                for ( int k = 0; k < num_ghost_cells; ++k )
                                {
                                    for ( int l = 0; l < this->nvar_; ++l )
                                    {
                                        std::vector<int> dof_indices;
                                        this->get_dofs_on_cell ( l, ghost_indices[k], dof_indices );
                                        MPI_Send ( vec2ptr ( dof_indices ), num_dofs_ghost[k * this->nvar_ + l], MPI_INT, i, 0, this->comm_ );
                                    }
                                }
                            }
                        }
                    }
                }
            }

            // Fix temporary Dof Ids on ghost layer to correct dof ids (this step might not always
            // be necessary but is essential to ensure correctness in special cases. See 
            // documentation file for more information

            max_dof_id = *( std::max_element ( this->numer_.begin ( ), this->numer_.end ( ) ) );

            std::vector<int> update_permutation ( max_dof_id + 1 );
            for ( size_t i = 0, e_i = update_permutation.size ( ); i != e_i; ++i )
                update_permutation[i] = i;

            for ( mesh::EntityIterator it = this->mesh_->begin ( this->tdim_ ),
                  e_it = this->mesh_->end ( this->tdim_ );
                  it != e_it;
                  ++it )
            {
                int subdomain_index;
                it->get ( "_sub_domain_", &subdomain_index );

                if ( subdomain_index != my_subdomain_ )
                {
                    int ghost_cell_index;
                    it->get ( "_remote_index_", &ghost_cell_index );

                    for ( size_t var = 0, e_var = this->nvar_; var != e_var; ++var )
                    {
                        // Get dofs from other subdomain
                        int size = this->fe_manager_->get_fe_on_cell ( it->index ( ), var )->get_nb_dof_on_cell ( );

                        std::vector<DofID> ghost_layer_dofs ( size );
                        for ( size_t i = 0, e_i = ghost_layer_dofs.size ( ); i != e_i; ++i )
                        {
                            ghost_layer_dofs[i] = numer_ghost[subdomain_index][ghost_cell_index][var][i];
                        }

                        // Get dofs on ghost layer from view of my subdomain
                        std::vector<DofID> critical_ghost_layer_dofs;
                        this->get_dofs_on_cell ( var, it->index ( ), critical_ghost_layer_dofs );

                        for ( size_t i = 0, e_i = critical_ghost_layer_dofs.size ( ); i != e_i; ++i )
                        {
                            const int cgld_i = critical_ghost_layer_dofs[i];
                            if ( cgld_i >= tmp_dof_offset + ndofs_on_sd_[my_subdomain_]
                                 || cgld_i < tmp_dof_offset )
                            {
                                assert ( cgld_i >= 0 && cgld_i < update_permutation.size ( ) );
                                update_permutation[cgld_i] = ghost_layer_dofs[i];
                            }
                        }
                    }
                }
            }

            this->apply_permutation ( update_permutation );

            LOG_DEBUG ( 2, " Third permutation done " );

            // Finaly calculate real dof_offset and correct numer_ field w.r.t. new offset

            my_dof_offset_ = 0;
            for ( size_t s = 0, e_s = my_subdomain_; s != e_s; ++s )
                my_dof_offset_ += ndofs_on_sd_[s];

            std::vector<int> old_dof_offsets ( number_of_subdomains_, 0 );
            for ( size_t s = 0, e_s = number_of_subdomains_; s != e_s; ++s )
                for ( size_t t = 0; t != s; ++t )
                    old_dof_offsets[s] += ndofs_with_ghost[t];

            std::vector<int> real_dof_offset ( number_of_subdomains_, 0 );
            for ( size_t s = 0, e_s = number_of_subdomains_; s != e_s; ++s )
                for ( size_t t = 0; t != s; ++t )
                    real_dof_offset[s] += ndofs_on_sd_[t];

            int size_of_ownerships = old_dof_offsets[number_of_subdomains_ - 1] + ndofs_with_ghost[number_of_subdomains_ - 1];
            std::vector<int> ownerships ( size_of_ownerships );

            for ( size_t s = 0, e_s = number_of_subdomains_ - 1; s != e_s; ++s )
                for ( size_t i = old_dof_offsets[s]; i < old_dof_offsets[s + 1]; ++i )
                    ownerships[i] = s;

            for ( size_t i = old_dof_offsets[number_of_subdomains_ - 1],
                  e_i = old_dof_offsets[number_of_subdomains_ - 1] + ndofs_with_ghost[number_of_subdomains_ - 1];
                  i < e_i;
                  ++i )
                ownerships[i] = number_of_subdomains_ - 1;

            max_dof_id = *( std::max_element ( this->numer_.begin ( ), this->numer_.end ( ) ) );

            std::vector<int> final_permutation ( max_dof_id + 1, -1 );
            for ( size_t i = 0, e_i = this->numer_.size ( ); i != e_i; ++i )
            {
                int owner = ownerships[this->numer_[i]];

                if ( owner != my_subdomain_ )
                    final_permutation[this->numer_[i]] = this->numer_[i] - old_dof_offsets[owner] + real_dof_offset[owner];
                else
                    final_permutation[this->numer_[i]] = this->numer_[i] - old_dof_offsets[my_subdomain_] + my_dof_offset_;
            }

            this->apply_permutation ( final_permutation );

            // Last check if Dof Ids are still greater than -1
            for ( size_t i = 0, e_i = this->numer_.size ( ); i != e_i; ++i )
            {
                assert ( this->numer_[i] >= 0 );
            }

            // Calculate number of dofs for each variable

            for ( size_t var = 0, e_var = this->nvar_; var != e_var; ++var )
            {
                int begin_offset = this->numer_offsets_cell_varloc_[var][0];
                int end_offset;

                if ( var + 1 < this->nvar_ )
                    end_offset = this->numer_offsets_cell_varloc_[var + 1][0];
                else
                    end_offset = this->numer_.size ( );

                std::vector<DofID> tmp ( end_offset - begin_offset );

                for ( size_t i = begin_offset; i < end_offset; ++i )
                    tmp[i - begin_offset] = this->numer_[i];

                std::sort ( tmp.begin ( ), tmp.end ( ) );

                std::vector<DofID>::iterator it = std::unique ( tmp.begin ( ), tmp.end ( ) );
                int tmp_size = it - tmp.begin ( );

                int hostile_dof = 0;
                for ( size_t i = 0, e_i = tmp_size; i != e_i; ++i )
                    if ( owner_of_dof ( tmp[i] ) != my_subdomain_ )
                        hostile_dof++;

                this->number_of_dofs_for_var_[var] -= hostile_dof;
            }

            // consecutive_numbering();
        }

        template<class DataType>
        void DofPartition<DataType>::consecutive_numbering ( )
        {
            // Fill vector local 2 global and global 2 local map
            local2global_.resize ( ndofs_incl_ghost_ );

            for ( size_t var = 0, e_var = this->nvar_; var != e_var; ++var )
            {
                int local_dof_cntr = 0;
                //First regular mesh without ghost layer
                for ( mesh::EntityIterator it = this->mesh_->begin ( this->tdim_ ),
                      e_it = this->mesh_->end ( this->tdim_ );
                      it != e_it;
                      ++it )
                {
                    int subdomain_index;
                    it->get ( "_sub_domain_", &subdomain_index );

                    if ( subdomain_index == my_subdomain_ )
                    {
                        std::vector<DofID> global_dofs;
                        this->get_dofs_on_cell ( var, it->index ( ), global_dofs );

                        for ( size_t i = 0, e_i = global_dofs.size ( ); i != e_i; ++i )
                        {
                            if ( global2local_.find ( global_dofs[i] ) == global2local_.end ( ) )
                            {
                                global2local_.insert ( std::pair<DofID, DofID>( global_dofs[i], local_dof_cntr ) );
                                ++local_dof_cntr;
                            }
                        }
                    }
                }
                // Next: Ghost layer
                for ( mesh::EntityIterator it = this->mesh_->begin ( this->tdim_ ),
                      e_it = this->mesh_->end ( this->tdim_ );
                      it != e_it;
                      ++it )
                {
                    int subdomain_index;
                    it->get ( "_sub_domain_", &subdomain_index );

                    if ( subdomain_index != my_subdomain_ )
                    {
                        std::vector<DofID> global_dofs;
                        this->get_dofs_on_cell ( var, it->index ( ), global_dofs );

                        for ( size_t i = 0, e_i = global_dofs.size ( ); i != e_i; ++i )
                        {
                            if ( global2local_.find ( global_dofs[i] ) == global2local_.end ( ) )
                            {
                                global2local_.insert ( std::pair<DofID, DofID>( global_dofs[i], local_dof_cntr ) );
                                ++local_dof_cntr;
                            }
                        }
                    }
                }
            }

            // Fill local2global
            std::map<DofID, DofID>::iterator it = global2local_.begin ( );
            while ( it != global2local_.end ( ) )
            {
                local2global_[( *it ).second] = ( *it ).first;
                ++it;
            }

            // Fill numer_sd_ field
            this->numer_sd_.resize ( this->numer_.size ( ), -1 );
            int offset = 0;
            for ( size_t var = 0, e_var = this->nvar_; var != e_var; ++var )
            {
                for ( mesh::EntityIterator it = this->mesh_->begin ( this->tdim_ ),
                      e_it = this->mesh_->end ( this->tdim_ );
                      it != e_it;
                      ++it )
                {
                    std::vector<DofID> global_dofs;
                    this->get_dofs_on_cell ( var, it->index ( ), global_dofs );

                    for ( size_t i = 0, e_i = global_dofs.size ( ); i != e_i; ++i )
                    {
                        DofID local_dof;
                        global2local ( global_dofs[i], &local_dof );
                        numer_sd_[i + offset] = local_dof;
                    }

                    offset += global_dofs.size ( );
                }
            }
        }

        template<class DataType>
        void DofPartition<DataType>::get_dofs_on_cell_sd ( int var, int cell_index, std::vector<DofID>& ids ) const
        {
            // finite element type

            FEType<DataType>& fe_type = *( this->fe_manager_->get_fe_on_cell ( cell_index, var ) );

            ids.resize ( fe_type.get_nb_dof_on_cell ( ) );

            // loop over DoFs

            for ( size_t i = 0, e_i = fe_type.get_nb_dof_on_cell ( ); i != e_i; ++i )
            {
                ids[i] = numer_sd_[this->numer_offsets_cell_varloc_[var][cell_index] + i];
            }
        }

        template<class DataType>
        void DofPartition<DataType>::apply_permutation_on_sd ( const std::vector<DofID>& permutation )
        {

            for ( size_t i = 0, e_i = numer_sd_.size ( ); i != e_i; ++i )
                numer_sd_[i] = permutation[numer_sd_[i]];

            // Fix local2global and global2local
            std::vector<DofID> l2g_backup ( local2global_ );
            for ( size_t i = 0, e_i = local2global_.size ( ); i != e_i; ++i )
                local2global_[i] = l2g_backup[permutation[local2global_[i]]];

            global2local_.clear ( );
            for ( size_t i = 0, e_i = local2global_.size ( ); i != e_i; ++i )
                global2local_.insert ( std::pair<DofID, DofID>( local2global_[i], i ) );
        }

        template<class DataType>
        void DofPartition<DataType>::GetPostProcessingStructure ( std::vector<DofID>& dof, std::vector<int>& dof_offset ) const
        {

            int nb_procs = -1;
            MPI_Comm_size ( comm_, &nb_procs );
            dof_offset.clear ( );
            dof_offset.resize ( nb_procs + 1, 0 );

            std::vector<DofID> dof_on_entity;
            SortedArray<DofID> dof_list;
            dof_list.reserve ( ndofs_on_sd ( my_subdomain_ ) );
            DofID id;

            for ( mesh::EntityIterator it = this->mesh_->begin ( this->tdim_ ),
                  e_it = this->mesh_->end ( this->tdim_ );
                  it != e_it;
                  ++it )
            {

                for ( size_t var = 0, e_var = this->nvar_; var != e_var; ++var )
                {
                    dof_on_entity.clear ( );
                    this->get_dofs_on_cell ( var, it->index ( ), dof_on_entity );
                    for ( size_t i = 0, e_i = dof_on_entity.size ( ); i != e_i; ++i )
                    {
                        id = dof_on_entity[i];
                        if ( !dof_list.find_insert ( id ) )
                        {
                            for ( size_t j = this->owner_of_dof ( id ) + 1; j < dof_offset.size ( ); ++j )
                            {
                                ++dof_offset[j];
                            }
                        }
                    }
                }
            }
            dof.clear ( );
            dof = dof_list.data ( );
        }

        // template instantiation
        template class DofPartition<double>;
        template class DofPartition<float>;

    }
} // namespace hiflow
