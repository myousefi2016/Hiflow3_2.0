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

#include "visualization.h"

#include <sstream>

#include "mesh/writer.h"
#include "common/macros.h"
#include "mesh/attributes.h"
#include "mesh/mesh_builder.h"
#include "mesh/mesh_db_view.h"

#include <mpi.h>

/// @author Martin Baumann, Thomas Gengenbach

namespace hiflow
{

    template<class DataType>
    Visualization<DataType>::Visualization ( )
    : visualize_all_attributes_ ( false )
    {
        LOG_INFO ( "WARNING", "This class is deprecated. Use CellVisualization if possible" );
        mesh_ = NULL;
        dof_ = NULL;
        fe_manager_ = NULL;
        filename_ = "visualization.vtu";
        var_names_.clear ( );
    }

    template<class DataType>
    Visualization<DataType>::~Visualization ( )
    {

    }

    template<class DataType>
    void Visualization<DataType>::set_var_names ( const std::vector<std::string>& names )
    {
        var_names_ = names;
        LOG_INFO ( "Variable names", string_from_range ( names.begin ( ), names.end ( ) ) );
        assert ( var_names_.size ( ) == fe_manager_->get_nb_var ( ) );
    }

    template<class DataType>
    void Visualization<DataType>::clear_var_names ( )
    {
        var_names_.clear ( );
    }

    template<class DataType>
    void Visualization<DataType>::set_visualize_all_attributes ( bool flag )
    {
        visualize_all_attributes_ = flag;
    }

    template<class DataType>
    void Visualization<DataType>::Visualize ( const std::vector<double>& visu_vec )
    {
        int num_partitions = -1;
        MPI_Comm_size ( dof_->get_mpi_comm ( ), &num_partitions );

        // checks
        interminable_assert ( mesh_ != NULL );
        interminable_assert ( dof_ != NULL );
        interminable_assert ( fe_manager_ != NULL );

        mesh::MeshPtr mesh_ptr;

        int size_of_vis_vector = 0;
        // size correct for sequential case must be adopted for parallel case
        size_of_vis_vector = mesh_->num_entities ( 0 ); // size of Q1 vector

        mesh::VtkWriter writer;
        mesh::PVtkWriter pwriter ( dof_->get_mpi_comm ( ) );

        if ( num_partitions > 1 )
        {
            mesh::TDim tdim = mesh_->tdim ( );

            if ( mesh_->has_attribute ( "_remote_index_", tdim ) )
            {
                // create new mesh without ghost layers
                mesh::MeshDbViewBuilder builder ( ( static_cast < mesh::MeshDbView* > ( mesh_ ) )->get_db ( ) );

                std::vector<int> index_map;

                // get __remote_indices__ attribute
                mesh::AttributePtr rem_indices_attr = mesh_->get_attribute ( "_remote_index_", tdim );

                // loop over cells that are not ghosts
                for ( mesh::EntityIterator cell_it = mesh_->begin ( tdim ); cell_it != mesh_->end ( tdim ); ++cell_it )
                {
                    if ( rem_indices_attr->get_int_value ( cell_it->index ( ) ) == -1 )
                    {
                        // loop over vertices of cell
                        std::vector<mesh::MeshBuilder::VertexHandle> vertex_handle;
                        for ( mesh::IncidentEntityIterator inc_vert_it = cell_it->begin_incident ( 0 ); inc_vert_it != cell_it->end_incident ( 0 ); ++inc_vert_it )
                        {
                            std::vector<double> coords;
                            inc_vert_it->get_coordinates ( coords );
                            vertex_handle.push_back ( builder.add_vertex ( coords ) );
                        }
                        builder.add_entity ( tdim, vertex_handle );
                        index_map.push_back ( cell_it->index ( ) );
                    }
                }
                mesh_ptr = builder.build ( );
                mesh::AttributePtr index_map_attr ( new mesh::IntAttribute ( index_map ) );
                mesh_ptr->add_attribute ( "__index_map__", mesh_->tdim ( ), index_map_attr );
            }
        }

        // Write all attributes in the mesh if the flag is set.
        if ( visualize_all_attributes_ )
        {
            if ( num_partitions == 1 )
            {
                writer.add_all_attributes ( *mesh_, true );
            }
            else
            {
                pwriter.add_all_attributes ( *mesh_, true );
            }
        }

        // parallel case
        if ( mesh_ptr.get ( ) != 0 )
        {
            Mesh* old_mesh = mesh_; // TODO: are we sure old_mesh has not been deleted here?
            // set mesh_ to new mesh without ghostcells in
            mesh_ = mesh_ptr.get ( );
            // adopt number of vertices to new mesh
            size_of_vis_vector = mesh_->num_entities ( 0 );

            // transfer cell attributes
            std::vector< std::string > cell_attr_names = old_mesh->get_attribute_names ( mesh_->tdim ( ) );

            mesh::AttributePtr index_map_attr = mesh_->get_attribute ( "__index_map__", mesh_->tdim ( ) );

            for ( std::vector< std::string >::const_iterator it = cell_attr_names.begin ( ),
                  end_it = cell_attr_names.end ( ); it != end_it; ++it )
            {
                mesh::AttributePtr mapped_attr (
                                                 new mesh::InheritedAttribute ( old_mesh->get_attribute ( *it, mesh_->tdim ( ) ),
                                                                                index_map_attr ) );
                mesh_->add_attribute ( *it, mesh_->tdim ( ), mapped_attr );
            }
        }
        // discontinuous ready visualization approach:
        // all cells may contribute to each vertex (sum) and afterward,
        // the sum is divided by number of contributions

        // visualization vectors

        // loop over variables
        for ( int var = 0; var < fe_manager_->get_nb_var ( ); ++var )
        {
            std::vector<double> temp_field ( size_of_vis_vector, 0.0 );
            std::vector<double> temp_weight ( size_of_vis_vector, 0.0 );

            // loop over cells of mesh without ghost cells
            MeshEntityIterator cell = mesh_->begin ( mesh_->tdim ( ) );

            while ( cell != mesh_->end ( mesh_->tdim ( ) ) )
            {
                std::vector<int> temp;
                // loop over vertices of cell
                MeshIncidentEntityIterator vertex_it = cell->begin_incident ( 0 );
                for ( int vertex = 0; vertex < cell->num_vertices ( ); ++vertex )
                {
                    // temp contains the dof id of vertex

                    // NB: This function is expected to get the dof id via the
                    // coordinates, if this is not the case anymore, the
                    // indeces of the newly created mesh have to be mapped
                    // to the old one.
                    // dof_ uses original mesh to find correct dof_ids
                    dof_->get_dofs_on_subentity ( var, cell->index ( ), 0, vertex, temp );
                    assert ( temp.size ( ) == 1 );

                    int dof_id = temp[0];
                    int vertex_index = vertex_it->index ( );

                    temp_field [vertex_index] += visu_vec[dof_id];
                    temp_weight[vertex_index]++;

                    ++vertex_it;
                }
                ++cell;
            }

            for ( int j = 0; j < temp_field.size ( ); ++j )
            {
                assert ( temp_weight.at ( j ) != 0.0 );
                temp_field.at ( j ) /= temp_weight.at ( j );
            }

            // write visualization file
            mesh::AttributePtr field_attribute ( new mesh::DoubleAttribute ( temp_field ) );
            std::string name;
            if ( var_names_.size ( ) == 0 )
            {
                std::stringstream var_str;
                var_str << var;
                name = std::string ( "val" ) + var_str.str ( );
            }
            else
            {
                name = var_names_.at ( var );
            }
            mesh_->add_attribute ( name.c_str ( ), 0, field_attribute );

            if ( num_partitions == 1 )
                writer.add_attribute ( name.c_str ( ), 0 );
            else
                pwriter.add_attribute ( name.c_str ( ), 0 );
        }

        if ( num_partitions == 1 )
        {
            writer.write ( filename_.c_str ( ), *mesh_ );
        }
        else
        {
            pwriter.write ( filename_.c_str ( ), *mesh_ );
        }
    }

    template class Visualization<double>;
    template class Visualization<float>;

} // namespace hiflow
