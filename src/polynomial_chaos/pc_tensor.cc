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

/// \author Michael Schick

#include "polynomial_chaos/pc_tensor.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <algorithm>
#include "common/macros.h"
#include "common/pointers.h"

namespace hiflow
{
    namespace polynomialchaos
    {

        PCTensor::PCTensor ( )
        {
            level_ = -1;
            q_ = -1;
            numproc_ = 0;
            my_rank_ = -1;
            comm_ = MPI_COMM_NULL;
            tensor_eps_ = 1.e-15;
        }

        PCTensor::~PCTensor ( )
        {
            int is_finalized;
            MPI_Finalized ( &is_finalized );
            if ( !is_finalized )
            {
                if ( this->comm_ != MPI_COMM_NULL )
                    MPI_Comm_free ( &this->comm_ );
            }
        }

        void PCTensor::SetMPIComm ( const MPI_Comm& comm )
        {
            assert ( comm != MPI_COMM_NULL );
            if ( this->comm_ != MPI_COMM_NULL )
                MPI_Comm_free ( &this->comm_ );
            int info = MPI_Comm_rank ( comm, &my_rank_ );
            assert ( info == MPI_SUCCESS );
            assert ( my_rank_ >= 0 );
            info = MPI_Comm_size ( comm, &numproc_ );
            assert ( info == MPI_SUCCESS );
            assert ( numproc_ > 0 );
            assert ( my_rank_ < numproc_ );

            info = MPI_Comm_split ( comm, 0, my_rank_, &comm_ );
            assert ( info == MPI_SUCCESS );
        }

        void PCTensor::ComputeOwnership ( PCBasis const& pcbasis )
        {
            ownership_.resize ( pcbasis.No ( ) + 1 );

            for ( int p = 0; p < pcbasis.No ( ) + 1; ++p )
            {
                for ( int k = 0; k < pc_size_[p]; ++k )
                    ownership_[p][k] = k % numproc_;
            }
        }

        void PCTensor::ComputeDistributedModesTable ( PCBasis const& pcbasis )
        {
            hostile_modes_needed_.resize ( pcbasis.No ( ) + 1 );
            hostile_modes_send_.resize ( pcbasis.No ( ) + 1 );
            noff_modes_.resize ( pcbasis.No ( ) + 1 );
            nmodes_local_.resize ( pcbasis.No ( ) + 1 );
            g2l_.resize ( pcbasis.No ( ) + 1 );
            l2g_.resize ( pcbasis.No ( ) + 1 );

            for ( int p = 0; p < pcbasis.No ( ) + 1; ++p )
                hostile_modes_needed_[p].resize ( numproc_ );

            for ( int p = 0; p < pcbasis.No ( ) + 1; ++p )
            {
                int mode_cntr = 0;
                for ( int i = 0; i < pc_size_[p]; ++i )
                {
                    if ( ownership_[p][i] == my_rank_ )
                    {
                        mode_cntr++;
                    }
                }
                nmodes_local_[p] = mode_cntr;

                std::vector<int> breakpoints = k_breakpoints_[p];
                std::vector<std::vector<int> > tmp;
                tmp.resize ( numproc_ );

                for ( int bp = 1; bp <= pc_size_[p]; ++bp )
                {
                    for ( int pos = breakpoints[bp - 1]; pos < breakpoints[bp]; ++pos )
                    {
                        std::vector<int> idx = glo2loc_[p][pos];
                        if ( ownership_[p][idx[2]] == my_rank_ && ownership_[p][idx[1]] != my_rank_ )
                        {
                            tmp[ownership_[p][idx[1]]].push_back ( idx[1] );
                        }
                    }
                }

                for ( int rank = 0; rank < numproc_; ++rank )
                {
                    hostile_modes_needed_[p][rank] = tmp[rank];
                    std::sort ( hostile_modes_needed_[p][rank].begin ( ), hostile_modes_needed_[p][rank].end ( ) );
                    std::vector<int>::iterator it;
                    it = std::unique ( hostile_modes_needed_[p][rank].begin ( ), hostile_modes_needed_[p][rank].end ( ) );
                    hostile_modes_needed_[p][rank].resize ( std::distance ( hostile_modes_needed_[p][rank].begin ( ), it ) );
                }

                int** sizes = new int* [numproc_];
                for ( int n = 0; n < numproc_; ++n )
                    sizes[n] = new int[numproc_];

                noff_modes_[p] = 0;
                for ( int n = 0; n < numproc_; ++n )
                {
                    sizes[my_rank_][n] = hostile_modes_needed_[p][n].size ( );
                    noff_modes_[p] += hostile_modes_needed_[p][n].size ( );
                }

                for ( int s = 0; s < numproc_; ++s )
                    MPI_Bcast ( sizes[s], numproc_, MPI_INT, s, comm_ );

                hostile_modes_send_[p].resize ( numproc_ );
                for ( int n = 0; n < numproc_; ++n )
                    hostile_modes_send_[p][n].resize ( sizes[n][my_rank_] );

                int nsend = 0;
                for ( int rank = 0; rank < numproc_; ++rank )
                {
                    if ( rank != my_rank_ && sizes[my_rank_][rank] > 0 )
                    {
                        nsend++;
                    }
                }

                std::vector<MPI_Request> request ( nsend );
                std::vector<MPI_Status> status ( nsend );

                int send_cntr = 0;
                for ( int rank = 0; rank < numproc_; ++rank )
                {
                    if ( rank != my_rank_ && sizes[my_rank_][rank] > 0 )
                    {
                        int tag = 99;
                        MPI_Isend ( &hostile_modes_needed_[p][rank][0], sizes[my_rank_][rank], MPI_INT, rank, tag, comm_, &request[send_cntr] );
                        send_cntr++;
                    }
                }

                for ( int rank = 0; rank < numproc_; ++rank )
                {
                    if ( rank != my_rank_ && sizes[rank][my_rank_] > 0 )
                    {
                        int tag = 99;
                        MPI_Status status;
                        MPI_Recv ( &hostile_modes_send_[p][rank][0], sizes[rank][my_rank_], MPI_INT, rank, tag, comm_, &status );
                    }
                }

                int info_send = MPI_Waitall ( nsend, &request.front ( ), &status.front ( ) );
                assert ( info_send == MPI_SUCCESS );

                for ( int n = 0; n < numproc_; ++n )
                    delete[] sizes[n];

                delete[] sizes;

                int cntr = 0;

                l2g_[p].resize ( nmodes_local_[p] + noff_modes_[p] );

                for ( int i = 0; i<static_cast < int > ( ownership_[p].size ( ) ); ++i )
                {
                    if ( ownership_[p][i] == my_rank_ )
                    {
                        g2l_[p][i] = cntr;
                        l2g_[p][cntr] = i;
                        cntr++;
                    }
                }

                for ( int rank = 0; rank < numproc_; ++rank )
                {
                    if ( rank != my_rank_ )
                    {
                        for ( int i = 0; i<static_cast < int > ( hostile_modes_needed_[p][rank].size ( ) ); ++i )
                        {
                            g2l_[p][hostile_modes_needed_[p][rank][i]] = cntr;
                            l2g_[p][cntr] = hostile_modes_needed_[p][rank][i];
                            cntr++;
                        }
                    }
                }
            } // for p=0,p<No...
        }

        bool pair_compare ( std::pair<int, int> const& lhs, std::pair<int, int> const& rhs )
        {
            return (lhs.second > rhs.second );
        }

        void PCTensor::ComputeTensor ( int input_degree, PCBasis const& pcbasis )
        {
            q_ = input_degree;
            int N = pcbasis.N ( );
            level_ = pcbasis.No ( );
            pc_size_.resize ( pcbasis.No ( ) + 1 );
            pos_linear_.resize ( pcbasis.No ( ) + 1 );
            k_breakpoints_.resize ( pcbasis.No ( ) + 1 );
            k_linear_breakpoints_.resize ( pcbasis.No ( ) + 1 );
            glo2loc_.resize ( pcbasis.No ( ) + 1 );
            val_.resize ( pcbasis.No ( ) + 1 );
            nnz_.resize ( pcbasis.No ( ) + 1 );

            std::vector<int> pc_size_input ( pcbasis.No ( ) + 1 );

            for ( int p = 0; p < pcbasis.No ( ) + 1; ++p )
            {
                k_breakpoints_[p].push_back ( 0 );
                k_linear_breakpoints_[p].push_back ( 0 );
                pc_size_[p] = pcbasis.SubDim ( p );
                pc_size_input[p] = pcbasis.SubDim ( p );
            }

            //int pc_size_q = 1;
            for ( int p = 0; p < pcbasis.No ( ) + 1; ++p )
            {
                if ( q_ < p ) pc_size_input[p] = pc_size_[q_];
            }

            square_norm_.resize ( pc_size_[pcbasis.No ( )] );

            for ( int k = 0; k < pcbasis.Dim ( ); ++k )
            {
                double square_norm = pcbasis.ComputeIntegral3rdOrder ( k, k, 0 );
                assert ( square_norm > 0.0 );
                square_norm_[k] = square_norm;
            }

            std::vector<int> tmp ( 3 );
            double val = 0.0;
            for ( int p = 0; p < pcbasis.No ( ) + 1; ++p )
            {
                for ( int k = 0; k < pc_size_[p]; ++k )
                {
                    for ( int j = 0; j < pc_size_[p]; ++j )
                    {
                        for ( int i = 0; i < pc_size_input[p]; ++i )
                        {
                            tmp[0] = i;
                            tmp[1] = j;
                            tmp[2] = k;
                            val = pcbasis.ComputeIntegral3rdOrder ( i, j, k );
                            if ( fabs ( val / square_norm_[k] ) > tensor_eps_ )
                            {
                                glo2loc_[p].push_back ( tmp );
                                val_[p].push_back ( val / square_norm_[k] );

                                if ( i <= N )
                                {
                                    pos_linear_[p].push_back ( glo2loc_[p].size ( ) - 1 );
                                }
                            }
                        }
                    }
                    k_breakpoints_[p].push_back ( val_[p].size ( ) );
                    k_linear_breakpoints_[p].push_back ( pos_linear_[p].size ( ) );
                }
                nnz_[p] = val_[p].size ( );
            }

            if ( comm_ != MPI_COMM_NULL )
            {
                ComputeOwnership ( pcbasis );
                ComputeDistributedModesTable ( pcbasis );
            }
        }

    }
}
