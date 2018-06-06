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

#include "numbering_strategy.h"
#include "degree_of_freedom.h"
#include "fem/femanager.h"
#include "mesh/entity.h"

/// \author Michael Schick<br>Martin Baumann

namespace hiflow
{
    namespace doffem
    {

        template<class DataType>
        void NumberingStrategy<DataType>::initialize ( DegreeOfFreedom<DataType>& dof, std::vector<DofID>& numer,
                                                       std::vector<std::vector<int> >& numer_offsets_cell_varloc,
                                                       DofInterpolation& dof_interpolation )
        {
            dof_ = &dof;
            mesh_ = &( dof_->get_mesh ( ) );
            fe_manager_ = &( dof_->get_fe_manager ( ) );
            numer_ = &numer;
            numer_offsets_cell_varloc_ = &numer_offsets_cell_varloc;
            dof_interpolation_ = &dof_interpolation;

            tdim_ = mesh_->tdim ( );
            nvar_ = fe_manager_->get_nb_var ( );
        }

        /// \param description is an optional parameter that should describe what
        ///                    the permutation represents

        template<class DataType>
        void NumberingStrategy<DataType>::apply_permutation ( const std::vector<DofID>& permutation,
                                                              const std::string& description )
        {
            // apply permutation to mapl2g

            // DoF IDs are used in numer_ only

            for ( size_t i = 0, e_i = numer_->size ( ); i != e_i; ++i )
                ( *numer_ )[i] = permutation[( *numer_ )[i]];

            // apply permutation to DofInterpolation

            dof_interpolation_->apply_permutation ( permutation );

            // calculate number of dofs, as this could have changed

            update_number_of_dofs ( description );
        }

        /// \param description is an optional parameter that should describe what
        ///                    the permutation represents

        template<class DataType>
        void NumberingStrategy<DataType>::update_number_of_dofs ( const std::string& description )
        {
            // Calculate number of DoFs

            number_of_dofs_total_ = 0;
            for ( size_t i = 0, e_i = numer_->size ( ); i != e_i; ++i )
                if ( ( *numer_ )[i] > number_of_dofs_total_ )
                    number_of_dofs_total_ = ( *numer_ )[i];
            ++number_of_dofs_total_;

            // Calculate number of Dofs for each variable

            number_of_dofs_for_var_.resize ( nvar_, 0 );
            for ( size_t var = 0; var != nvar_; ++var )
            {
                int begin_offset = ( *numer_offsets_cell_varloc_ )[var][0];
                int end_offset;

                if ( var < nvar_ - 1 )
                    end_offset = ( *numer_offsets_cell_varloc_ )[var + 1][0];
                else
                    end_offset = numer_->size ( );

                for ( size_t i = begin_offset; i < end_offset; ++i )
                    if ( ( *numer_ )[i] > number_of_dofs_for_var_[var] )
                        number_of_dofs_for_var_[var] = ( *numer_ )[i];

                for ( size_t j = 0; j != var; ++j )
                    number_of_dofs_for_var_[var] -= number_of_dofs_for_var_[j];

                ++number_of_dofs_for_var_[var];
            }

#if 0
            int rank = -1;
            MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
            if ( ( number_of_dofs_total_old != number_of_dofs_total_ ) &&
                 ( rank == 0 ) )
                std::cout << "#Dofs: " << number_of_dofs_total_
                    << " " << description << std::endl;
#endif
        }

        template<class DataType>
        void NumberingStrategy<DataType>::get_points ( const mesh::Entity& entity, std::vector<Coord>& points )
        {
            points.reserve ( points.size ( ) + entity.num_vertices ( ) );
            for ( size_t p = 0; p != entity.num_vertices ( ); ++p )
            {
                std::vector<DataType> coords;
                entity.get_coordinates ( coords );
                points.push_back ( coords );
            }
        }

        // template instantiation
        template class NumberingStrategy<double>;
        template class NumberingStrategy<float>;

    }
}
