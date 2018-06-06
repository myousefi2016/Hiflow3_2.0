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

#include "periodic_boundary_conditions_cartesian.h"
#include "common/macros.h"
#include "dof/degree_of_freedom.h"
#include "mesh/iterator.h"
#include "mesh/mesh.h"
#include "mesh/types.h"
#include "space/dirichlet_boundary_conditions.h"

/// @author Martin Baumann

namespace hiflow
{

    template<class DataType>
    void PeriodicBoundaryConditionsCartesian<DataType>::compute_conditions ( VectorSpace<DataType>& space )
    {

        int rank = -1;
        MPI_Comm_rank ( space.get_mpi_comm ( ), &rank );

        const mesh::Mesh& mesh = space.mesh ( );
        const mesh::TDim tdim = mesh.tdim ( );
        const mesh::GDim gdim = mesh.gdim ( );

        // *****************************************************
        // 1. Initalization

        x_coord_.clear ( );
        y_coord_.clear ( );
        z_coord_.clear ( );

        this->corresponding_dof_.clear ( );

        // *****************************************************
        // 2. Collect DoFs on boundaries

        // extract boundary of mesh
        mesh::MeshPtr boundary_mesh = mesh.extract_boundary_mesh ( );

        // TODO: should be possible to use own communicator here
        int size = -1;
        MPI_Comm_size ( space.get_mpi_comm ( ), &size );
        const bool is_sequential = ( size == 1 );
        if ( !is_sequential )
            assert ( mesh.has_attribute ( "_sub_domain_", tdim ) );

        // loop over boundary faces
        for ( mesh::EntityIterator it_boundary = boundary_mesh->begin ( tdim - 1 );
              it_boundary != boundary_mesh->end ( tdim - 1 );
              ++it_boundary )
        {
            // get id of boundary face
            const mesh::Id boundary_id = it_boundary->id ( );

            // get index number where to store boundary entity data
            int face_number;
            const bool check = mesh.find_entity ( tdim - 1, boundary_id, &face_number );
            assert ( check );

            // Get the face to be able to access to the data associated with the face
            mesh::Entity face = mesh.get_entity ( tdim - 1, face_number );

            // reset the cell because it was changed in the assert above
            mesh::IncidentEntityIterator cell = face.begin_incident ( tdim );

            // get owner of cell
            int subdomain;
            if ( !is_sequential )
                cell->get ( "_sub_domain_", &subdomain );
            else
                subdomain = rank;

            if ( subdomain == rank )
            {
                // get material number on face
                typename PeriodicBoundaryConditionsCartesian<DataType>::MaterialNumber mat_num = face.get_material_number ( );

                // get dofs and coords on face for a given variable
                std::vector<doffem::DofID> dofs_on_face;
                std::vector<std::vector<DataType> > coords_on_face;

                for ( int var = 0; var < space.get_nb_var ( ); ++var )
                {
                    // Get the dofs indices and coords on the face
                    find_dofs_and_coords_on_face<DataType>( it_boundary,
                            space,
                            var,
                            dofs_on_face,
                            coords_on_face );

                    // store information for later construction of corresponding DoFs
                    for ( int i = 0; i < dofs_on_face.size ( ); ++i )
                    {
                        if ( is_periodic_boundary ( mat_num ) )
                        {
                            bool new_entry; // true if entry hasn't been considered
                            new_entry = insert_dof_id_entry ( mat_num, var, dofs_on_face[i] );

                            if ( new_entry == true )
                            {
                                x_coord_.insert ( std::make_pair ( dofs_on_face[i], coords_on_face[i][0] ) );
                                y_coord_.insert ( std::make_pair ( dofs_on_face[i], coords_on_face[i][1] ) );
                                if ( gdim == 3 )
                                    z_coord_.insert ( std::make_pair ( dofs_on_face[i], coords_on_face[i][2] ) );
                            }
                        } // if (is_periodic_boundary(mat_num))
                    } // for (int i=0; i<dofs_on_face.size(); ++i)
                } // for (int var=0; var<space.get_nb_var(); ++var)
            } // if (subdomain == rank)
        } // for (EntityIterator it_boundary = boundary_mesh->begin(tdim-1);

        // *****************************************************
        // 3. Detect corresponding DoFs

        int verbatim = 0; // 0 -> no output
        // 1 -> with output

        if ( verbatim > 0 )
        {
            std::cout << "====================================================================" << std::endl;
            std::cout << std::endl;
            std::cout << "Periodic Boundaries: Matching DOFs" << std::endl;
        }

        // Iterate over (mat, mat) tuples

        typename std::map<typename PeriodicBoundaryConditionsCartesian<DataType>::MaterialNumber, typename PeriodicBoundaryConditionsCartesian<DataType>::MaterialNumber>::const_iterator mat_tuple_it = this->boundary_descriptor_.begin ( );
        while ( mat_tuple_it != this->boundary_descriptor_.end ( ) )
        {

            typename PeriodicBoundaryConditionsCartesian<DataType>::MaterialNumber mat_a = mat_tuple_it->first;
            typename PeriodicBoundaryConditionsCartesian<DataType>::MaterialNumber mat_b = mat_tuple_it->second;

            // not to handle each material tuple twice ...

            if ( mat_a < mat_b )
            {

                // Iterate over vars

                //for (int var=0; var<space.get_nb_var()-1; ++var) //no pressure DOFs for Navier-Stokes
                for ( int var = 0; var < space.get_nb_var ( ); ++var )
                {

                    if ( verbatim > 0 )
                    {
                        std::cout << std::endl;
                        std::cout << "MAT_A " << mat_a << "\t MAT_B " << mat_b << "\t var " << var << std::endl;
                    }

                    // Vectors containing DoF Ids on mat_a and mat_b

                    typename std::map < std::pair<typename PeriodicBoundaryConditionsCartesian<DataType>::MaterialNumber, int>, std::vector<doffem::DofID> >::const_iterator it;

                    it = this->dof_id_map_.find ( std::make_pair ( mat_a, var ) );
                    interminable_assert ( it != this->dof_id_map_.end ( ) );
                    std::vector<doffem::DofID> const& dofs_a = it->second;

                    it = this->dof_id_map_.find ( std::make_pair ( mat_b, var ) );
                    interminable_assert ( it != this->dof_id_map_.end ( ) );
                    std::vector<doffem::DofID> const& dofs_b = it->second;

                    interminable_assert ( dofs_a.size ( ) > 1 );
                    interminable_assert ( dofs_a.size ( ) == dofs_b.size ( ) );

                    // Determine match mode, i.e. ordinate

                    // Principle:
                    //  x ordinates constant => (y,z) used for DoF match
                    //  y ordinates constant => (x,z) used for DoF match
                    //  z ordinates constant => (x,y) used for DoF match

                    std::map<doffem::DofID, double>* ordinate_1 = 0;
                    std::map<doffem::DofID, double>* ordinate_2 = 0;

                    if ( z_coord_.size ( ) == 0 )
                    {
                        // 2D case
                        if ( ( x_coord_[dofs_a[0]] == x_coord_[dofs_a[1]] ) &&
                             ( y_coord_[dofs_a[0]] != y_coord_[dofs_a[1]] ) )
                            ordinate_1 = &y_coord_;
                        else if ( ( x_coord_[dofs_a[0]] != x_coord_[dofs_a[1]] ) &&
                                  ( y_coord_[dofs_a[0]] == y_coord_[dofs_a[1]] ) )
                            ordinate_1 = &x_coord_;
                        else
                            quit_program ( );
                    }
                    else
                    {
                        // 3D case
                        interminable_assert ( z_coord_.size ( ) > 1 );
                        interminable_assert ( x_coord_.size ( ) == y_coord_.size ( ) );
                        interminable_assert ( y_coord_.size ( ) == z_coord_.size ( ) );

                        bool x_changing = false;
                        bool y_changing = false;
                        bool z_changing = false;

                        unsigned number_changings = 0;

                        for ( int i = 1; i < dofs_a.size ( ); ++i )
                        {
                            if ( ( x_coord_[dofs_a[i]] != x_coord_[dofs_a[0]] ) &&
                                 ( x_changing == false ) )
                            {
                                x_changing = true;
                                ++number_changings;
                            }
                            if ( ( y_coord_[dofs_a[i]] != y_coord_[dofs_a[0]] ) &&
                                 ( y_changing == false ) )
                            {
                                y_changing = true;
                                ++number_changings;
                            }
                            if ( ( z_coord_[dofs_a[i]] != z_coord_[dofs_a[0]] ) &&
                                 ( z_changing == false ) )
                            {
                                z_changing = true;
                                ++number_changings;
                            }
                            if ( number_changings >= 2 )
                                break;
                        }

                        interminable_assert ( number_changings == 2 );

                        if ( x_changing && y_changing )
                        {
                            ordinate_1 = &x_coord_;
                            ordinate_2 = &y_coord_;
                        }
                        else if ( y_changing && z_changing )
                        {
                            ordinate_1 = &y_coord_;
                            ordinate_2 = &z_coord_;
                        }
                        else if ( x_changing && z_changing )
                        {
                            ordinate_1 = &x_coord_;
                            ordinate_2 = &z_coord_;
                        }
                    }

                    // Matching process

                    for ( int i = 0; i < dofs_a.size ( ); ++i )
                    {

                        std::vector<doffem::DofID> dof_vector;
                        dof_vector.reserve ( dofs_a.size ( ) * dofs_b.size ( ) );

                        if ( z_coord_.size ( ) == 0 )
                        {
                            // 2D case
                            double ord = ordinate_1->find ( dofs_a[i] )->second;

                            // find j such that ordinate of dofs_b[j] = ord
                            for ( int j = 0; j < dofs_b.size ( ); ++j )
                            {
                                if ( ordinate_1->find ( dofs_b[j] )->second == ord )
                                {
                                    dof_vector.push_back ( dofs_b[j] );
                                    break;
                                }
                            }
                        }
                        else
                        {
                            // 3D case
                            for ( int i = 0; i < dofs_a.size ( ); ++i )
                            {
                                double ord_1 = ordinate_1->find ( dofs_a[i] )->second;
                                double ord_2 = ordinate_2->find ( dofs_a[i] )->second;

                                // find j such that ordinate of dofs_b[j] = ord
                                for ( int j = 0; j < dofs_b.size ( ); ++j )
                                {
                                    if ( ordinate_1->find ( dofs_b[j] )->second == ord_1 )
                                    {
                                        if ( ordinate_2->find ( dofs_b[j] )->second == ord_2 )
                                        {
                                            dof_vector.push_back ( dofs_b[j] );
                                            break;
                                        }
                                    }
                                }
                            }
                        } // if (z_coord_.size() == 0)

                        // insert tuple
                        interminable_assert ( dof_vector.size ( ) >= 0 );
                        this->corresponding_dof_.insert ( std::make_pair ( dofs_a[i], dof_vector ) );

                        if ( verbatim > 0 )
                        {
                            std::cout << dofs_a[i] << ": \t x = " << x_coord_[dofs_a[i]] << "\t"
                                    "y = " << y_coord_[dofs_a[i]] << "\t"
                                    << dof_vector.at ( 0 ) << ": \t x = " << x_coord_[dof_vector.at ( 0 )] << "\t"
                                    "y = " << y_coord_[dof_vector.at ( 0 )] << std::endl;
                        } // if (verbatim > 0)

                    } // if (z_coord_.size() == 0)
                } // for (int var=0; ...
            } // if (mat_a < mat_b)

            ++mat_tuple_it;

        } // while (mat_tuple_it != _pbh->end())

        this->handle_multiple_correspondences ( );

        if ( verbatim > 0 )
        {
            std::cout << std::endl;
            std::cout << "====================================================================" << std::endl;
        }
    }

} // namespace hiflow
