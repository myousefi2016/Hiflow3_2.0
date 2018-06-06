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

/// \author Chandramowli Subramanian, Simon Gawlok

#ifndef HIFLOW_INTEGRATION_H_
#    define HIFLOW_INTEGRATION_H_

#    include <cmath>
#    include <cstdio>
#    include <string>

#    include "common/sorted_array.h"
#    include "mesh/mesh.h"
#    include "linear_algebra/seq_dense_matrix.h"
#    include "quadrature/quadrature.h"
#    include "space/solution.h"
#    include "space/vector_space.h"
#    include "common/csv_writer.h"

/// If defined, sparsity structure is written to structure.csv. Only meaningful 
/// in sequential mode.
//#define OUTPUT_STRUCT

namespace hiflow
{

    /// Set up non zero structure in order to assemble matrix.
    /// (row/col) pairs correspond to components in the matrix to set
    /// Checks for duplicates.
    /// @param space vector space
    /// @param rows_diag global row numbers for diagonal block
    /// @param cols_diag global column numbers for diagonal block
    /// @param rows_offdiag global row numbers for offdiagonal block
    /// @param cols_offdiag global column numbers for offdiagonal block
    /// \param[in] coupling_vars 2D array indicating the coupling of vars.
    ///                          Rows (first index) belong to test variables,
    ///                          columns (second index) belong to trial variables.
    ///                          If entry (i, j) is set to true, test variable i
    ///                          couples with trial variable j, and the corresponding
    ///                          block is contained in the sparsity structure. Otherwise,
    ///                          the block is skipped, resulting in a sparser structure.
    ///                          If this argument is not passed or is empty, full coupling
    ///                          of all variables is assumed. All rows and columns need
    ///                          to have the size of space.get_nb_var().

    template<class DataType>
    inline void InitStructure ( const VectorSpace<DataType>& space,
                                std::vector<int>* rows_diag,
                                std::vector<int>* cols_diag,
                                std::vector<int>* rows_offdiag,
                                std::vector<int>* cols_offdiag,
                                std::vector<std::vector<bool> > *coupling_vars )
    {

        // Assert correct size of coupling_vars

        assert ( coupling_vars->size ( ) == space.get_nb_var ( ) );
        for ( int i = 0, i_e = space.get_nb_var ( ); i != i_e; ++i )
        {
            assert ( ( *coupling_vars )[i].size ( ) == space.get_nb_var ( ) );
        }

#    ifdef OUTPUT_STRUCT
        // CSV writer for visual check of sparsity structure
        CSVWriter<int> writer ( "sparsity.csv" );

        std::vector<std::string> names;
        names.push_back ( "i" );
        names.push_back ( "j" );
        names.push_back ( "val" );

        writer.Init ( names );

        std::vector<int> values ( 3, -1 );
#    endif

        // total number of own dofs
        int ndof_total = space.dof ( ).ndofs_on_sd ( space.dof ( ).get_my_subdomain ( ) );

        // a set of columns for every row
        // std::vector< std::set<int> > raw_struct_diag(ndof_total);
        // std::vector< std::set<int> > raw_struct_offdiag(ndof_total);
        std::vector< SortedArray<int> > raw_struct_diag ( ndof_total );
        std::vector< SortedArray<int> > raw_struct_offdiag ( ndof_total );

        std::vector<int> dof_ind_test, dof_ind_trial;
        int local_dof_i;

        // loop over every cell (including ghost cells)
        typename VectorSpace<DataType>::MeshEntityIterator mesh_it = space.mesh ( ).begin ( space.get_dim ( ) );
        typename VectorSpace<DataType>::MeshEntityIterator e_mesh_it = space.mesh ( ).end ( space.get_dim ( ) );
        while ( mesh_it != e_mesh_it )
        {
            // loop over test variables
            for ( int test_var = 0, tv_e = space.get_nb_var ( ); test_var != tv_e; ++test_var )
            {
                // get dof indices for test variable
                space.GetDofIndices ( test_var, *mesh_it, &dof_ind_test );

                // loop over trial variables
                for ( int trial_var = 0, vt_e = space.get_nb_var ( ); trial_var != vt_e; ++trial_var )
                {

                    // check whether test_var and trial_var couple
                    if ( ( *coupling_vars )[test_var][trial_var] )
                    {

                        // get dof indices for trial variable
                        space.GetDofIndices ( trial_var, *mesh_it, &dof_ind_trial );

                        // detect couplings
                        for ( size_t i = 0, i_e = dof_ind_test.size ( ); i != i_e; ++i )
                        {
                            const int di_i = dof_ind_test[i];

                            // if my row
                            if ( space.dof ( ).is_dof_on_sd ( di_i ) )
                            {

                                space.dof ( ).global2local ( di_i, &local_dof_i );

                                for ( size_t j = 0, j_e = dof_ind_trial.size ( ); j != j_e; ++j )
                                {
                                    const int di_j = dof_ind_trial[j];

                                    // diagonal coupling (my col)
                                    if ( space.dof ( ).is_dof_on_sd ( di_j ) )
                                    {
                                        raw_struct_diag[local_dof_i].find_insert ( di_j );
                                    }
                                    else
                                    {
                                        // nondiagonal coupling (not my col)
                                        raw_struct_offdiag[local_dof_i].find_insert ( di_j );
                                    }
                                } // endif my row

                            } // for (int j=0;...
                        } // for (int i=0;...
                    }
                }
            }
            // next cell
            ++mesh_it;
        } // while (mesh_it != ...

#    ifdef OUTPUT_STRUCT
        for ( size_t k = 0, k_e = raw_struct_diag.size ( ); k != k_e; ++k )
        {
            values[0] = k;
            for ( size_t l = 0, l_e = raw_struct_diag[k].size ( ); l != l_e; ++l )
            {
                values[1] = raw_struct_diag[k][l];
                values[2] = 1;
                writer.write ( values );
            }
        }

        // compute nnz for nondiagonal block
        for ( size_t k = 0, k_e = raw_struct_offdiag.size ( ); k != k_e; ++k )
        {
            values[0] = k;
            for ( size_t l = 0, l_e = raw_struct_offdiag[k].size ( ); l != l_e; ++l )
            {
                values[1] = raw_struct_offdiag[k][l];
                values[2] = 1;
                writer.write ( values );
            }
        }
#    endif

        // compute nnz for diagonal block
        int nnz_diag = 0;
        for ( size_t k = 0, k_e = raw_struct_diag.size ( ); k != k_e; ++k )
        {
            nnz_diag += raw_struct_diag[k].size ( );
        }

        // compute nnz for nondiagonal block
        int nnz_offdiag = 0;
        for ( size_t k = 0, k_e = raw_struct_offdiag.size ( ); k != k_e; ++k )
        {
            nnz_offdiag += raw_struct_offdiag[k].size ( );
        }

        // now create rows, cols
        rows_diag->resize ( nnz_diag );
        cols_diag->resize ( nnz_diag );
        rows_offdiag->resize ( nnz_offdiag );
        cols_offdiag->resize ( nnz_offdiag );

        // now create row/col indices for diagonal block
        int ind = 0;
        int global_dof_i = 0;

        for ( size_t i = 0; i != ndof_total; ++i )
        {
            // now going through row i
            space.dof ( ).local2global ( i, &global_dof_i );

            for ( SortedArray<int>::const_iterator
                  it = raw_struct_diag[i].begin ( ),
                  end = raw_struct_diag[i].end ( );
                  it != end; ++it )
            {

                ( *rows_diag )[ind] = global_dof_i;
                ( *cols_diag )[ind] = *it;

                // increase index
                ++ind;
            }
        }
        assert ( ind == nnz_diag );

        // now create row/col indices for offdiagonal block
        ind = 0;
        for ( size_t i = 0; i != ndof_total; ++i )
        {
            // now going through row i
            space.dof ( ).local2global ( i, &global_dof_i );

            for ( SortedArray<int>::const_iterator
                  it = raw_struct_offdiag[i].begin ( ),
                  end = raw_struct_offdiag[i].end ( );
                  it != end; ++it )
            {

                ( *rows_offdiag )[ind] = global_dof_i;
                ( *cols_offdiag )[ind] = *it;
                // increase index
                ++ind;
            }
        }
        assert ( ind == nnz_offdiag );
    }

}

#endif  // HIFLOW_INTEGRATION_H_
