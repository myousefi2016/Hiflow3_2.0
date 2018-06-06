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

#include "point_evaluator.h"

namespace hiflow
{

    template<class DataType>
    PointEvaluator<DataType>::PointEvaluator ( const VectorSpace<DataType>& space )
    : space_ ( space )
    {
        trial_cells_.resize ( 0 );
    }

    template<class DataType>
    bool PointEvaluator<DataType>::evaluate_fun ( const EvalFunction& fun,
                                                  const std::vector<Coordinate>& point,
                                                  DataType& value ) const
    {
        std::vector<int> cell_index;
        return evaluate_fun ( fun, point, value, cell_index );
    }

    template<class DataType>
    bool PointEvaluator<DataType>::evaluate_fun ( const EvalFunction& fun,
                                                  const std::vector<Coordinate>& point,
                                                  DataType& value,
                                                  std::vector<int>& cell_index ) const
    {
        value = 0.0;
        mesh::MeshPtr mesh = space_.meshPtr ( );
        mesh::GeometricSearchPtr search = mesh->get_geometric_search ( );

        //      std::vector<int> cell_index;
        if ( trial_cells_.size ( ) > 0 )
        {
            search->find_cell ( point, trial_cells_, cell_index );
        }
        else
        {
            search->find_cell ( point, cell_index );
        }

        int value_count = 0;
        if ( !cell_index.empty ( ) )
        {
            for ( int i = 0; i < cell_index.size ( ); ++i )
            {
                Entity cell = mesh->get_entity ( mesh->tdim ( ), cell_index[i] );
                int rem_ind = -100;
                cell.get<int>( "_remote_index_", &rem_ind );
                if ( rem_ind == -1 )
                {
                    const typename doffem::CellTransformation<DataType>& ct = space_.GetCellTransformation ( cell );
                    std::vector<Coordinate> ref_coords;
                    bool found = ct.contains_physical_point ( point, &ref_coords );
                    if ( found )
                    {
                        std::vector<DataType> temp_value ( 1, 1.e32 );
                        fun ( cell, ref_coords, temp_value );
                        value += temp_value[0];
                        ++value_count;
                    }
                }
            }

            if ( value_count > 0 )
            {
                value /= static_cast < DataType > ( value_count );
                return true;
            }
            return false;
        }
        else
        {
            return false;
        }
    }

    template<class DataType>
    bool PointEvaluator<DataType>::evaluate_fun_global ( const EvalFunction& fun,
                                                         const std::vector<Coordinate>& point,
                                                         DataType& value,
                                                         const MPI_Comm& comm ) const
    {
        std::vector<int> cell_index;
        return evaluate_fun_global ( fun, point, value, cell_index, comm );
    }

    template<class DataType>
    bool PointEvaluator<DataType>::evaluate_fun_global ( const EvalFunction& fun,
                                                         const std::vector<Coordinate>& point,
                                                         DataType& value,
                                                         std::vector<int>& cell_index,
                                                         const MPI_Comm& comm ) const
    {
        int rank;
        MPI_Comm_rank ( comm, &rank );

        DataType has_point;
        has_point = static_cast < DataType > ( evaluate_fun ( fun, point, value, cell_index ) );

        DataType sum_data[2];

        sum_data[0] = 0.;
        sum_data[1] = 0.;

        DataType origin_data[2];
        // write recv/send data in array for less communication
        origin_data[0] = value;
        origin_data[1] = has_point;
        MPI_Allreduce ( origin_data,
                        sum_data,
                        2,
                        mpi_data_type<DataType>::get_type ( ),
                        MPI_SUM,
                        comm );
        if ( sum_data[1] > 0.0 )
        {
            origin_data[0] = sum_data[0] / sum_data[1];
            origin_data[1] = 1.;
        }
        else
        {
            // no process found the point
            origin_data[0] = 0.;
            origin_data[1] = 0.;
        }
        value = origin_data[0];
        return static_cast < bool > ( origin_data[1] );
    }

    template<class DataType>
    void PointEvaluator<DataType>::set_trial_cells ( const std::vector<int>& trial_cells )
    {
        trial_cells_.resize ( trial_cells.size ( ), 0 );
        for ( int l = 0; l < trial_cells.size ( ); ++l )
        {
            trial_cells_[l] = trial_cells[l];
        }
    }

    // TODO: As geometric quantities are given fixed in double precision,
    // compatibility problems occur using float here.
    //  template class PointEvaluator<float>;
    template class PointEvaluator<double>;
}
