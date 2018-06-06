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

#include "fe_interpolation.h"
#include "common/log.h"
#include "common/vector_algebra.h"
#include "dof/degree_of_freedom.h"
#include "fem/fetype.h"
#include "mesh/geometric_search.h"
#include "assembly/assembly_assistant_values.h"
#include "assembly/assembly.h"

#include <cmath>
#include <utility>
#include <string>
#include <vector>
#include <mpi.h>
#include <sstream>
#include <algorithm>

const int DEBUG_LEVEL = 1;

namespace hiflow
{
    ////////////////////////////////////////////////////
    /////// SpacePatchInterpolation /////////////////////
    ////////////////////////////////////////////////////

    template<class LAD, int DIM>
    FEInterpolation<LAD, DIM>::FEInterpolation ( )
    : eps_ ( 1e-12 ),
    num_var_ ( 0 ),
    initialized_ ( false ),
    from_space_ ( NULL ),
    to_space_ ( NULL ),
    tdim_ ( 0 ),
    gdim_ ( 0 ),
    rank_ ( -1 )
    {
        this->degrees_.resize ( num_var_ );
        this->is_cg_.resize ( num_var_ );
    }

    template<class LAD, int DIM>
    void FEInterpolation<LAD, DIM>::init ( VectorSpace<DataType>* from_space, VectorSpace<DataType> * to_space )
    {
        this->clear ( );
        this->from_space_ = from_space;
        this->to_space_ = to_space;
        this->from_mesh_ = from_space->meshPtr ( );
        this->to_mesh_ = to_space->meshPtr ( );
        this->tdim_ = from_mesh_->tdim ( );
        this->gdim_ = to_mesh_->gdim ( );
        this->num_var_ = from_space->get_nb_var ( );

        assert ( this->from_mesh_ != NULL );

        const MPI_Comm& comm = from_space->get_mpi_comm ( );

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
            this->is_cg_[v] = from_space->fe_manager ( ).get_ca ( v );
            this->degrees_[v] = from_space->fe_manager ( ).get_fe_on_cell ( 0, v )->get_fe_deg ( ) * 2;
        }
    }

    template<class LAD, int DIM>
    void FEInterpolation<LAD, DIM>::clear ( )
    {
        this->num_var_ = 0;
        this->eps_ = 1e-12;
        this->degrees_.clear ( );
        this->is_cg_.clear ( );
        this->tdim_ = 0;
        this->gdim_ = 0;
        this->from_mesh_ = NULL;
        this->from_space_ = NULL;
        this->to_mesh_ = NULL;
        this->to_space_ = NULL;
        this->rank_ = 0;
        this->initialized_ = false;
    }
    /*
    template class FEInterpolation<LADescriptorCoupledD, 1>;
    template class FEInterpolation<LADescriptorCoupledD, 2>;
    template class FEInterpolation<LADescriptorCoupledD, 3>;
#ifdef WITH_HYPRE
    template class FEInterpolation<LADescriptorHypreD, 1>;
    template class FEInterpolation<LADescriptorHypreD, 2>;
    template class FEInterpolation<LADescriptorHypreD, 3>;
#endif
     */
    ///////////////////////////////////////////////////////
    //////// FEInterpolationLagrange //////////////////////
    ///////////////////////////////////////////////////////

    template<class LAD, int DIM>
    FEInterpolationLagrange<LAD, DIM>::FEInterpolationLagrange ( )
    : FEInterpolation<LAD, DIM>( )
    {
    }

    template<class LAD, int DIM>
    void FEInterpolationLagrange<LAD, DIM>::clear ( )
    {
        FEInterpolation<LAD, DIM>::clear ( );
        this->dof_weights_.clear ( );
    }

    template<class LAD, int DIM>
    void FEInterpolationLagrange<LAD, DIM>::init ( VectorSpace<DataType>* from_space, VectorSpace<DataType>* to_space )
    {
        FEInterpolation<LAD, DIM>::init ( from_space, to_space );
        this->create_dof_mapping ( );
        this->initialized_ = true;
    }

    template<class LAD, int DIM>
    bool FEInterpolationLagrange<LAD, DIM>::check_space ( ) const
    {
        // TODO

        return true;
    }

    template<class LAD, int DIM>
    void FEInterpolationLagrange<LAD, DIM>::create_dof_mapping ( )
    {
        this->dof_weights_.resize ( this->num_var_ );

        // setup search object for from-mesh
        GridGeometricSearch grid_search ( this->from_mesh_ );

        // loop over all variables
        for ( int v = 0; v<this->num_var_; ++v )
        {
            int num_dofs_to = this->to_space_->dof ( ).get_nb_dofs ( v );

            // loop over all cells in to-mesh
            for ( mesh::EntityIterator it = this->to_mesh_->begin ( this->tdim_ ),
                  end_it = this->to_mesh_->end ( this->tdim_ ); it != end_it; ++it )
            {
                std::vector<int> global_dof_ids_on_cell_to;
                this->to_space_->GetDofIndices ( v, *it, &global_dof_ids_on_cell_to );
                int num_dofs_on_cell_to = global_dof_ids_on_cell_to.size ( );
                std::vector<DataType> values;
                values.resize ( num_dofs_on_cell_to, 0. );

                std::vector< Coord > coords_to;
                this->to_space_->dof ( ).get_coord_on_cell ( v, it->index ( ), coords_to );

                std::vector<int> trial_cells_from;

                // loop over to-dofs on current cell
                for ( int j = 0; j < num_dofs_on_cell_to; ++j )
                {
                    int global_j = global_dof_ids_on_cell_to[j];

                    // consider only locally owned dofs
                    if ( this->rank_ != this->to_space_->dof ( ).owner_of_dof ( global_j ) )
                    {
                        continue;
                    }

                    // check if dof has been already considered
                    if ( this->dof_weights_[v].find ( global_j ) != this->dof_weights_[v].end ( ) )
                    {
                        continue;
                    }

                    // insert empty map
                    std::map<int, DataType> map_j;

                    // search for cells in from_mesh which contain to-dof
                    std::vector<int> cells_from;
                    grid_search.find_cell ( coords_to[j], trial_cells_from, cells_from );

                    // if no cell is found, continue with next to-dof
                    if ( cells_from.size ( ) == 0 )
                    {
                        continue;
                    }

                    // augment trial cells by found cells
                    trial_cells_from.insert ( trial_cells_from.end ( ), cells_from.begin ( ), cells_from.end ( ) );

                    //if (print) std::cout << "dof_j " << global_j << " with coords " << coords_to[j][0] << ", " << coords_to[j][1] << std::endl;

                    // loop over all shape functions of from-space with support in cells_j
                    int num_cells_from = cells_from.size ( );
                    if ( num_cells_from > 1 )
                    {
                        num_cells_from = 1;
                    }

                    for ( int c = 0; c < num_cells_from; ++c )
                    {
                        // get from-cell entity and fe-type on cell
                        int cell_index_from = cells_from[c];
                        Entity from_cell = this->from_mesh_->get_entity ( this->tdim_, cell_index_from );
                        std::vector< const FEType<DataType>* > fe_type_from;
                        fe_type_from.push_back ( this->from_space_->fe_manager ( ).get_fe_on_cell ( cell_index_from, v ) );

                        // setup evalshapefunctions object
                        EvalShapeFunctions<DIM, DataType> shape_fun ( fe_type_from );

                        // get reference coordinates for to-dof in current from-cell
                        CellTransformation<DataType>* cell_trafo = this->from_space_->fe_manager ( ).get_cell_transformation ( cell_index_from );
                        Vec<DIM, DataType> ref_coord;
                        if ( this->tdim_ == 1 )
                        {
                            cell_trafo->inverse ( coords_to[j][0], ref_coord[0] );
                        }
                        else if ( this->tdim_ == 2 )
                        {
                            cell_trafo->inverse ( coords_to[j][0], coords_to[j][1], ref_coord[0], ref_coord[1] );
                            //if (print) std::cout << "ref coord " << ref_coord[0] << " " << ref_coord[1] << std::endl;
                        }
                        else if ( this->tdim_ == 3 )
                        {
                            cell_trafo->inverse ( coords_to[j][0], coords_to[j][1], coords_to[j][2], ref_coord[0], ref_coord[1], ref_coord[2] );
                        }

                        // loop over all shape functions on from-cell
                        std::vector<int> global_dof_ids_on_cell_from;
                        this->from_space_->GetDofIndices ( v, from_cell, &global_dof_ids_on_cell_from );
                        int num_dofs_on_cell_from = global_dof_ids_on_cell_from.size ( );

                        // evaluate from-shape_fun's at to-dof
                        std::vector<DataType> shape_vals;
                        shape_fun ( -1, ref_coord, shape_vals );

                        assert ( shape_vals.size ( ) == num_dofs_on_cell_from );

                        for ( int i = 0; i < num_dofs_on_cell_from; ++i )
                        {
                            int global_i = global_dof_ids_on_cell_from[i];

                            // store result in interpolator-mapping
                            if ( std::abs ( shape_vals[i] ) > this->eps_ )
                            {
                                //if (print) std::cout << "interpolating dof " << global_i << " with shape value " << shape_vals[i] << std::endl;
                                map_j.insert ( std::pair<int, DataType>( global_i, shape_vals[i] ) );
                            }
                        }
                    }
                    this->dof_weights_[v].insert ( std::pair<int, std::map<int, DataType> > ( global_j, map_j ) );
                }
            }
        }
    }

    template<class LAD, int DIM>
    void FEInterpolationLagrange<LAD, DIM>::interpolate ( const VectorType& from_vec, VectorType& to_vec ) const
    {
        if ( !this->initialized_ )
        {
            std::cout << "FEInterpolation is not initialized ! " << std::endl;
            exit ( -1 );
        }

        to_vec.Zeros ( );
        int in_begin = from_vec.ownership_begin ( );
        int in_end = from_vec.ownership_end ( );
        int out_begin = to_vec.ownership_begin ( );
        int out_end = to_vec.ownership_end ( );

        /*
        std::cout << "  Ownership range in-vector:    " << in_begin << " - " << in_end << std::endl;
        std::cout << "  Ownership range out-vector:    " << out_begin << " - " << out_end << std::endl;
         */

        // loop over all variables
        for ( int v = 0; v<this->num_var_; ++v )
        {
            // loop over all dofs for current variable
            for ( std::map<int, std::map<int, double> >::const_iterator it_j = this->dof_weights_[v].begin ( ); it_j != this->dof_weights_[v].end ( ); ++it_j )
            {
                int j = it_j->first;
                if ( this->rank_ != this->to_space_->dof ( ).owner_of_dof ( j ) )
                {
                    continue;
                }

                // loop over interpolating dofs of from-space
                for ( std::map<int, double >::const_iterator it_i = it_j->second.begin ( ); it_i != it_j->second.end ( ); ++it_i )
                {
                    std::vector<int> i ( 1 );
                    i[0] = it_i->first;
                    DataType coeff = it_i->second;

                    std::vector<DataType> in_val ( 1, 0. );

                    LOG_DEBUG ( 2, "dof j = " << j << " interpolated by dof i = " << i[0] << " with weight " << coeff );

                    from_vec.GetValues ( &i[0], 1, &in_val[0] );
                    to_vec.Add ( j, coeff * in_val[0] );
                }
            }
        }
        to_vec.Update ( );
        interpolate_constrained_vector ( *this->to_space_, to_vec );
        to_vec.Update ( );
    }

    template class FEInterpolationLagrange<LADescriptorCoupledD, 1>;
    template class FEInterpolationLagrange<LADescriptorCoupledD, 2>;
    template class FEInterpolationLagrange<LADescriptorCoupledD, 3>;
#ifdef WITH_HYPRE
    template class FEInterpolationLagrange<LADescriptorHypreD, 1>;
    template class FEInterpolationLagrange<LADescriptorHypreD, 2>;
    template class FEInterpolationLagrange<LADescriptorHypreD, 3>;
#endif

}
