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

///
/// \file eval_cyl_fe.cc
/// \brief Function for evaluating a FE function in cylindrical coordinates
///
/// \author Martin Baumann, Philipp Gerstner
///

#include "eval_cyl_fe.h"

void EvalCylFeFunction::operator() ( const Entity& cell,
        const std::vector<double>& ref_coords,
        std::vector<double>& values ) const
{

    const int gdim = space_.mesh ( ).gdim ( );
    const int num_points = ref_coords.size ( ) / gdim;

    std::vector<int> global_dof_ids_phi;
    std::vector<int> global_dof_ids_r;
    std::vector<int> global_dof_ids_z;

    space_.GetDofIndices ( 0, cell, &global_dof_ids_phi );
    space_.GetDofIndices ( 1, cell, &global_dof_ids_r );
    if ( DIM == 3 )
        space_.GetDofIndices ( 2, cell, &global_dof_ids_z );

    const int num_dofs = global_dof_ids_phi.size ( );

    interminable_assert ( num_dofs == global_dof_ids_r.size ( ) );

    std::vector< std::vector<double> > shape_fun ( num_points, std::vector<double>( num_dofs, 1.e13 ) );

    std::vector<double> pt ( gdim, 0. );

    int k = 0;
    for ( int i = 0; i < num_points; ++i )
    {
        for ( int c = 0; c < gdim; ++c )
        {
            pt[c] = ref_coords[k++];
        }
        space_.fe_manager ( ).get_fe_on_cell ( cell.index ( ), var_ )->N ( pt, shape_fun[i] );
    }

    std::vector< doffem::DegreeOfFreedom<DATATYPE>::Coord > coords;
    space_.dof ( ).get_coord_on_cell ( var_, cell.index ( ), coords );

    values.clear ( );
    values.resize ( num_points*DIM, 0. );

    for ( int j = 0; j < num_dofs; j++ )
    {
        {
            double dof_values_uphi;
            double dof_values_ur;
            double dof_values_uz;
            double dof_values_x;
            double dof_values_y;
            double dof_values_z;

            sol_.GetValues ( &global_dof_ids_phi[j], 1, &dof_values_uphi );
            sol_.GetValues ( &global_dof_ids_r[j], 1, &dof_values_ur );
            if ( DIM == 3 )
                sol_.GetValues ( &global_dof_ids_z[j], 1, &dof_values_uz );

            double phi = coords[j][0];
            dof_values_x = cos ( phi ) * dof_values_ur - sin ( phi ) * dof_values_uphi;
            dof_values_y = sin ( phi ) * dof_values_ur + cos ( phi ) * dof_values_uphi;
            if ( DIM == 3 )
                dof_values_z = dof_values_uz;

            for ( int i = 0; i < num_points; ++i )
            {
                values[DIM * i] += shape_fun.at ( i ).at ( j ) * dof_values_x;
                values[1 + DIM * i] += shape_fun.at ( i ).at ( j ) * dof_values_y;
                if ( DIM == 3 )
                    values[2 + DIM * i] += shape_fun.at ( i ).at ( j ) * dof_values_z;
            }
        }
    }
}
