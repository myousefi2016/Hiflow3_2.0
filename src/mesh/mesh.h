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

/// \author Thomas Gengenbach and Staffan Ronnas

#ifndef HIFLOW_MESH_MESH_H
#    define HIFLOW_MESH_MESH_H

#    include <map>
#    include <vector>
#    include <string>
#    include "mpi.h"
#    include "mesh/attributes.h"
#    include "mesh/types.h"
#    include "common/hdf5_tools.h"
#    include "mesh/interface.h"
#    include "mesh/periodicity_tools.h"

namespace MPI
{
    class Comm;
}

namespace hiflow
{
    namespace mesh
    {
        // Forward declarations
        class Connectivity;
        class Entity;
        class EntityIterator;
        class IncidentEntityIterator;
        class MeshDatabase;
        class RefinementTree;
        class Mesh;

        // Type of implementation

        enum IMPL
        {
            IMPL_DBVIEW,
            IMPL_REFINED,
            IMPL_P4EST
        };

        ///
        /// \brief Abstract base class for mesh implementations.
        ///

        class Mesh
        {
          public:
            /// \brief Constructor
            ///
            /// \param tdim  Maximal topological dimension of the entities in the mesh
            /// \param gdim  Geometrical dimension of the vertices in the mesh
            /// \param period Describing periodicity. Default: No periodicity
            inline Mesh ( TDim tdim, GDim gdim, std::vector<MasterSlave> period = std::vector<MasterSlave>( 0 ) );

            /// \brief Destructor
            inline virtual ~Mesh ( );

            /// \brief Access the topological dimension of the mesh
            ///
            /// \return  The topological dimension
            inline TDim tdim ( ) const;

            /// \brief Access the geometrical dimension of the mesh
            ///
            /// \return  The geometrical dimension
            inline GDim gdim ( ) const;

            /// \brief Add a new attribute to the mesh
            ///
            /// \param name  Name of attribute
            /// \param tdim  Topological dimension of entities for which attribute is defined
            /// \param attribute A pointer to the attribute
            inline void add_attribute ( const std::string& name, TDim tdim, AttributePtr attribute );

            /// \brief Access an attribute by name
            ///
            /// \param name  Name of attribute
            /// \param tdim  Topological dimension of entities for which attribute is defined
            /// \return A pointer to the attribute
            inline AttributePtr get_attribute ( const std::string& name, TDim tdim ) const;

            /// \brief Get all attribute names
            ///
            /// \param tdim  Topological dimension of entities for which
            /// attribute names are required
            /// \return Vector with attribute names as strings
            inline std::vector<std::string> get_attribute_names ( TDim Tdim ) const;

            /// \brief Access the value of an attribute
            ///
            /// \param name  Name of attribute
            /// \param tdim  Topological dimension of entities for which attribute is defined
            /// \param index Index of entity
            /// \param[out] value  Pointer where value should be stored
            template<typename T>
            void get_attribute_value ( const std::string& name, TDim tdim, EntityNumber index, T* value ) const;

            /// \brief Set the value of an attribute
            ///
            /// \param name  Name of attribute
            /// \param tdim  Topological dimension of entities for which attribute is defined
            /// \param index Index of entity
            /// \param value New value of attribute for the entity
            template<typename T>
            void set_attribute_value ( const std::string& name, TDim tdim, EntityNumber index, const T& value ) const;

            /// \brief Check whether a given attribute exists in the Mesh
            ///
            /// \param name  Name of attribute
            /// \param tdim  Topological dimension of entities for which attribute is defined
            /// \return true if an attribute with the given name and dimension exists in the Mesh
            inline bool has_attribute ( const std::string& name, TDim tdim ) const;

            /// \brief Obtain iterator to the entities of the mesh
            /// topological dimension in the mesh
            ///
            /// \param entity_dim  The topological dimension of the entities to iterate over
            /// \return  An iterator to the first entity of dimension entity_dim in the mesh
            virtual EntityIterator begin ( TDim entity_dim ) const = 0;

            /// \brief Obtain iterator to the entities of the mesh
            /// topological dimension in the mesh
            ///
            /// \param entity_dim  The topological dimension of the entities to iterate over
            /// \return  An iterator to one-past-the-last entity of dimension entity_dim in the mesh
            virtual EntityIterator end ( TDim entity_dim ) const = 0;

            /// \brief Obtain iterator to the first entity of dimension entity_dim which is incident to a given entity.
            virtual IncidentEntityIterator begin_incident ( const Entity& entity, TDim entity_dim ) const = 0;

            /// \brief Obtain iterator to one-past-the-last entity of dimension entity_dim which is incident to a given entity.
            virtual IncidentEntityIterator end_incident ( const Entity& entity, TDim entity_dim ) const = 0;

            /// \brief Access the id of an entity.
            virtual Id get_id ( TDim entity_dim, EntityNumber index ) const = 0;

            /// \brief Obtain the vertex id:s of an entity.
            virtual std::vector<Id> get_vertex_ids ( TDim entity_dim, EntityNumber index ) const = 0;

            /// \brief Access the material number of an entity
            virtual MaterialNumber get_material_number ( TDim entity_dim, EntityNumber index ) const = 0;

            /// \brief Change the material number of an entity
            virtual void set_material_number ( TDim entity_dim, EntityNumber index, MaterialNumber material ) = 0;

            /// \brief Access the coordinates of an entity.
            virtual std::vector<Coordinate> get_coordinates ( TDim entity_dim, EntityNumber index ) const = 0;

            /// \brief Access the id of the parent cell of a cell entity.
            virtual Id get_parent_cell_id ( EntityNumber cell_index ) const = 0;

            /// \brief Access the parent cell entity of a refined cell entity.
            /// \pre cell_has_parent(cell_index) == true
            virtual Entity get_parent_cell ( EntityNumber cell_index ) const = 0;

            /// \brief Check whether a cell has a parent cell.
            virtual bool cell_has_parent ( EntityNumber cell_index ) const = 0;

            /// \brief Obtain the id:s of the children of a cell entity.
            virtual std::vector<Id> get_children_cell_ids ( EntityNumber cell_index ) const = 0;

            virtual std::vector<EntityNumber> get_sibling_cell_indices ( EntityNumber cell_index ) const = 0;

            /// \brief Access the number of entities in the mesh
            ///
            /// \param entity_dim  The topological dimension of the entities
            /// \return  The number of entities of topological dimension entity_dim
            virtual EntityCount num_entities ( TDim entity_dim ) const = 0;

            /// \brief Access the number of entities incident to a given entity
            ///
            /// \param entity  An entity in the mesh
            /// \param entity_dim  The topological dimension of the entities
            /// \return  The number of entities of topological dimension entity_dim incident to the given entity
            virtual EntityCount num_incident_entities ( const Entity& entity, TDim entity_dim ) const = 0;

            /// \brief Obtain a view of a single entity in the mesh
            ///
            /// \param entity_dim  The topological dimension of the desired entity
            /// \param entity_number  The index of the entity in the mesh
            /// \return  An Entity object allowing access to the data associated with the entity
            virtual Entity get_entity ( TDim entity_dim, EntityNumber entity_number ) const = 0;

            /// \brief Look up an entity in the mesh by its id
            ///
            /// \param entity_dim  The topological dimension of the desired entity
            /// \param id  The id of the entity
            /// \param entity_number  An optional pointer to a location where the entity number should be stored if it is found.
            /// \return  true if the mesh contains an entity with the given dimension and id
            virtual bool find_entity ( TDim entity_dim, Id id, EntityNumber* entity_number ) const = 0;

            /// \brief Look up a vertex by its coordinates
            ///
            /// \param      coords  Coordinate vector of vertex
            /// \param[out] v_index Index of vertex if found
            /// \return true if vertex was found
            virtual bool find_vertex ( const std::vector<Coordinate>& coords, EntityNumber* v_index ) const = 0;

            /// \brief Refine the mesh uniformly
            ///
            /// Refines all the cells of the mesh with the default type of refinement (refinement_type = 0).
            /// \return  A shared pointer to a Mesh object containing the refined cells.
            virtual MeshPtr refine ( ) const = 0;

            /// \brief Refine the mesh with the provided refinement types.
            ///
            /// \param     refinements  refinement type of each cell, -1 for coarsening, other negative number for no action.
            /// \return  A shared pointer to a Mesh object containing the refined cells.
            virtual MeshPtr refine ( std::vector<int>& refinements ) const = 0;

            /// \brief Refine the mesh with given refinement types.
            ///
            /// Refines cells based on the value of a cell attribute.  The
            /// attribute should be of type int, and contain the desired
            /// refinement_type for each cell. -1 is used to indicate that
            /// no refinement should be performed on a cell.
            /// \param   attribute_name  Name of attribute containing the desired refinement type for each cell.
            /// \return  A shared pointer to a Mesh object containing the refined cells.
            virtual MeshPtr refine ( const std::string& attribute_name ) const = 0;

            virtual MeshPtr refine_uniform_seq ( int num_ref_steps ) const = 0;
            virtual MeshPtr refine_uniform_seq ( ) const = 0;

            /// \brief Replace a vertex to another location
            ///
            /// Replaces a vertex of the mesh to a desired destination.
            /// The user has to take care of the effects if further data
            /// structures have already been defined based on the mesh.
            /// For example the cell-transformations of the dof-module
            /// have to be reinitialized.
            /// \param  destination    New location of the vertex.
            /// \param  vertex_index   Index of the vertex to be moved.
            virtual void replace_vertex ( const std::vector<Coordinate>& destination, Id vertex_index ) = 0;

            /// \brief Move a vertex along a displacement vector
            ///
            /// Moves a vertex of the mesh along a displacement vector.
            /// The user has to take care of the effects if further data
            /// structures have already been defined based on the mesh.
            /// For example the cell-transformations of the dof-module
            /// have to be reinitialized.
            /// \param  displacement   Displacement of the vertex.
            /// \param  vertex_index   Index of the vertex to be displaced.
            virtual void move_vertex ( const std::vector<Coordinate>& displacement, Id vertex_index ) = 0;

            /// \brief Displaces all vertices of the mesh
            ///
            /// Moves all vertices of the mesh along respective
            /// displacement vectors.
            /// The user has to take care of the effects if further data
            /// structures have already been defined based on the mesh.
            /// For example the cell-transformations of the dof-module
            /// have to be reinitialized.
            /// \param  displacements    A vector of the size gdim * number_of_vertices.
            ///                          It contains the displacement of each vertex.
            virtual void move_vertices ( const std::vector<Coordinate>& displacements ) = 0;

            /// \brief Extracts the boundary of the mesh
            ///
            /// Identifies the facets (entities of dimension this->tdim()-1)
            /// that lie on the boundary, and creates a new mesh with these entities.
            /// \return A new Mesh object of topological dimension tdim()-1 containing the boundary facets.
            /// Ownership of the object is transferred to the caller.
            virtual MeshPtr extract_boundary_mesh ( ) const = 0;

            /// \brief Determines if facet is on the boundary
            /// \param[in] facet_index The index of the current facet
            /// \return True if it is a boundary facet, false otherwise
            virtual bool is_boundary_facet ( EntityNumber facet_index ) const = 0;

            /// \brief Get a geometric search for the mesh
            /// \return Smartpointer to the search instance
            virtual GeometricSearchPtr get_geometric_search ( ) = 0;

            /// \brief Reset the geometric search for the mesh
            virtual void reset_geometric_search ( ) = 0;

            virtual void set_refined_geometry_function ( RefinedGeometryFunction f ) = 0;
            virtual void set_recursive_refined_geometry_function ( RecursiveRefinedGeometryFunction f ) = 0;

            /// \brief compute global number of cells in mesh, shared by all processes in given communicator
            /// @param[in] comm MPI communicator
            /// @return number of cells
            virtual int num_global_cells ( const MPI_Comm& comm ) const = 0;

            /// \brief compute number of locally owned (i.e. no ghost) cells in mesh
            /// @return number of cells
            virtual int num_local_cells ( ) const = 0;

            /// \brief compute number of ghost cells in mesh
            /// @return number of cells
            virtual int num_ghost_cells ( ) const = 0;

            /// \brief return true if cell is locally ownder, i.e. not contained in ghost layer
            virtual bool cell_is_local ( EntityNumber index ) const = 0;

            /// \brief Save current state of mesh in file
            inline virtual void save ( std::string filename, const MPI_Comm& comm ) const;

            /// \brief Load current state of mesh in file
            inline virtual void load ( std::string filename, const MPI_Comm& comm );

            /// \brief copy mesh object from given reference
            inline virtual void copy_from ( const MeshPtr mesh );

            /// \brief copy data from reference mesh, including underlying data structures, to which the mesh has a pointer
            inline virtual void deep_copy_from ( const MeshPtr mesh );

            /// \brief Get periodicity
            /// \return Vector containing the information about the periodicity
            /// of the mesh
            inline virtual std::vector<MasterSlave> get_period ( ) const;

            /// \brief Check if mesh is periodic
            /// \return true if periodic else false
            inline virtual bool is_periodic ( ) const;

            inline virtual int get_ghost_layer_width ( ) const;

            inline virtual bool is_uniformly_coarsenable ( ) const;

          protected:

            mutable int ref_count_;
            ScopedArray<AttributeTable>::Type attributes_;
            std::vector<MasterSlave> period_;
            mutable int ghost_layer_width_;

          private:
            const TDim tdim_;
            const GDim gdim_;

            friend void intrusive_ptr_add_ref ( const Mesh* );
            friend void intrusive_ptr_release ( const Mesh* );
        };

        //////////////// Mesh Implementation ////////////////

        Mesh::Mesh ( TDim tdim, GDim gdim, std::vector<MasterSlave> period )
        : tdim_ ( tdim ), gdim_ ( gdim ), ref_count_ ( 0 ), attributes_ ( new AttributeTable[tdim + 1] ), period_ ( period ), ghost_layer_width_ ( 1 )
        {
        }

        Mesh::~Mesh ( )
        {
        }

        TDim Mesh::tdim ( ) const
        {
            return tdim_;
        }

        GDim Mesh::gdim ( ) const
        {
            return gdim_;
        }

        void Mesh::add_attribute ( const std::string& name, TDim tdim, AttributePtr attribute )
        {
            assert ( tdim >= 0 );
            assert ( tdim <= tdim_ );
            attributes_[tdim].add_attribute ( name, attribute );
        }

        AttributePtr Mesh::get_attribute ( const std::string& name, TDim tdim ) const
        {
            assert ( tdim >= 0 );
            assert ( tdim <= tdim_ );
            return attributes_[tdim].get_attribute ( name );
        }

        std::vector<std::string> Mesh::get_attribute_names ( TDim tdim ) const
        {
            assert ( tdim >= 0 );
            assert ( tdim <= tdim_ );
            std::vector<std::string> names ( attributes_[tdim].get_attribute_names ( ) );
            return names;
        }

        template<typename T>
        void Mesh::get_attribute_value ( const std::string& name, TDim tdim, EntityNumber index, T* value ) const
        {
            assert ( tdim >= 0 );
            assert ( tdim <= tdim_ );
            return attributes_[tdim].get ( name, index, value );
        }

        template<typename T>
        void Mesh::set_attribute_value ( const std::string& name, TDim tdim, EntityNumber index, const T& value ) const
        {
            assert ( tdim >= 0 );
            assert ( tdim <= tdim_ );
            return attributes_[tdim].set ( name, index, value );
        }

        bool Mesh::has_attribute ( const std::string& name, TDim tdim ) const
        {
            return attributes_[tdim].has_attribute ( name );
        }

        std::vector<MasterSlave> Mesh::get_period ( ) const
        {
            return period_;
        }

        bool Mesh::is_periodic ( ) const
        {
            return (period_.size ( ) > 0 );
        }

        int Mesh::get_ghost_layer_width ( ) const
        {
            return this->ghost_layer_width_;
        }

        bool Mesh::is_uniformly_coarsenable ( ) const
        {
            return false;
        }

        void Mesh::copy_from ( const MeshPtr mesh )
        {
            assert ( this->tdim ( ) == mesh->tdim ( ) );
            assert ( this->gdim ( ) == mesh->gdim ( ) );

            this->ghost_layer_width_ = mesh->get_ghost_layer_width ( );

            this->attributes_.reset ( new AttributeTable[mesh->tdim ( ) + 1] );

            // copy attribute table
            // Loop over all dimensions
            for ( int dim = 0; dim < mesh->tdim ( ) + 1; ++dim )
            {
                std::vector<std::string> names ( mesh->attributes_[dim].get_attribute_names ( ) );

                // Loop over all attributes
                for ( std::vector<std::string>::iterator it = names.begin ( ); it != names.end ( ); ++it )
                {
                    AttributePtr attr_ptr = mesh->attributes_[dim].get_attribute ( *it );

                    Attribute* attr = attr_ptr.get ( );
                    IntAttribute* int_attr;
                    DoubleAttribute* double_attr;
                    BoolAttribute* bool_attr;

                    int i = 0;
                    if ( ( int_attr = dynamic_cast < IntAttribute* > ( attr ) ) != 0 )
                    {
                        std::vector<int> int_att ( this->num_entities ( dim ) );
                        for ( std::vector<int>::iterator it_att = int_att.begin ( ); it_att < int_att.end ( ); ++it_att, ++i )
                        {
                            *it_att = int_attr->get_int_value ( i );
                        }

                        AttributePtr new_attribute = AttributePtr ( new IntAttribute ( int_att ) );
                        this->add_attribute ( *it, dim, new_attribute );

                    }
                    else if ( ( double_attr = dynamic_cast < DoubleAttribute* > ( attr ) ) != 0 )
                    {
                        std::vector<double> double_att ( this->num_entities ( dim ) );
                        for ( std::vector<double>::iterator it_att = double_att.begin ( ); it_att < double_att.end ( ); ++it_att, ++i )
                        {
                            *it_att = double_attr->get_double_value ( i );
                        }

                        AttributePtr new_attribute = AttributePtr ( new DoubleAttribute ( double_att ) );
                        this->add_attribute ( *it, dim, new_attribute );
                    }
                    else if ( ( bool_attr = dynamic_cast < BoolAttribute* > ( attr ) ) != 0 )
                    {
                        std::vector<bool> bool_att ( this->num_entities ( dim ) );
                        for ( std::vector<bool>::iterator it_att = bool_att.begin ( ); it_att < bool_att.end ( ); ++it_att, ++i )
                        {
                            *it_att = bool_attr->get_bool_value ( i );
                        }

                        AttributePtr new_attribute = AttributePtr ( new BoolAttribute ( bool_att ) );
                        this->add_attribute ( *it, dim, new_attribute );
                    }
                }
            }

            // copy period
            this->period_ = mesh->get_period ( );
        }

        void Mesh::deep_copy_from ( const MeshPtr mesh )
        {
            this->copy_from ( mesh );
        }

        void Mesh::save ( std::string filename, const MPI_Comm& comm ) const
        {
#    ifdef WITH_HDF5
            H5FilePtr file_ptr ( new H5File ( filename, "w", comm ) );

            //SETTING UP HDF5 GROUP
            std::stringstream groupname;
            groupname << "Mesh";
            H5GroupPtr group_ptr ( new H5Group ( file_ptr, groupname.str ( ), "w" ) );

            for ( int dim = 0; dim < tdim ( ) + 1; ++dim )
            {
                std::vector<std::string> names ( attributes_[dim].get_attribute_names ( ) );

                for ( std::vector<std::string>::iterator it = names.begin ( ); it != names.end ( ); ++it )
                {
                    AttributePtr attr_ptr = attributes_[dim].get_attribute ( *it );
                    Attribute* attr = attr_ptr.get ( );
                    IntAttribute* int_attr;
                    DoubleAttribute* double_attr;
                    BoolAttribute* bool_attr;

                    std::string att_prefix = "";
                    int i = 0;
                    if ( ( int_attr = dynamic_cast < IntAttribute* > ( attr ) ) != 0 )
                    {
                        att_prefix = "i";
                        std::vector<int> int_att ( this->num_entities ( dim ) );
                        for ( std::vector<int>::iterator it_att = int_att.begin ( ); it_att < int_att.end ( ); ++it_att, ++i )
                        {
                            *it_att = int_attr->get_int_value ( i );
                        }
                        write_array_parallel<int>( group_ptr, *it, int_att, comm );
                    }
                    else if ( ( double_attr = dynamic_cast < DoubleAttribute* > ( attr ) ) != 0 )
                    {
                        att_prefix = "d";
                        std::vector<double> double_att ( this->num_entities ( dim ) );
                        for ( std::vector<double>::iterator it_att = double_att.begin ( ); it_att < double_att.end ( ); ++it_att, ++i )
                        {
                            *it_att = double_attr->get_double_value ( i );
                        }
                        write_array_parallel<double>( group_ptr, *it, double_att, comm );
                    }
                    else if ( ( bool_attr = dynamic_cast < BoolAttribute* > ( attr ) ) != 0 )
                    {
                        att_prefix = "b";
                        std::vector<int> bool_att ( this->num_entities ( dim ) );
                        for ( std::vector<int>::iterator it_att = bool_att.begin ( ); it_att < bool_att.end ( ); ++it_att, ++i )
                        {
                            *it_att = bool_attr->get_bool_value ( i );
                        }
                        write_array_parallel<int>( group_ptr, *it, bool_att, comm );
                    }
                    std::stringstream new_att_name;
                    new_att_name << att_prefix << *it;
                    *it = new_att_name.str ( );
                }
                std::stringstream attr_names;
                attr_names << "attribute_names";
                if ( dim < 10 )
                    attr_names << "0";
                attr_names << dim;
                if ( names.size ( ) == 0 )
                    names.push_back ( "__empty__" );
                write_array_parallel ( group_ptr, attr_names.str ( ), names, comm );
            }

            // Write period
            //convert period to double vector
            std::vector<double> conv_period ( 4 * period_.size ( ) );
            for ( int k = 0; k < period_.size ( ); ++k )
            {
                conv_period[4 * k + 0] = period_[k].master ( );
                conv_period[4 * k + 1] = period_[k].slave ( );
                conv_period[4 * k + 2] = period_[k].h ( );
                conv_period[4 * k + 3] = ( double ) period_[k].index ( );
            }
            if ( period_.size ( ) == 0 ) // can't write empty vectors
                conv_period.push_back ( 0. );
            write_array_parallel<double>( group_ptr, "period", conv_period, comm );

#    else
            throw "HiFlow was not compiled with HDF5 support!\n";
#    endif
        }

        void Mesh::load ( std::string filename, const MPI_Comm& comm )
        {
#    ifdef WITH_HDF5
            H5FilePtr file_ptr ( new H5File ( filename, "r", comm ) );

            //SETTING UP HDF5 GROUP
            std::stringstream groupname;
            groupname << "Mesh";
            H5GroupPtr group_ptr ( new H5Group ( file_ptr, groupname.str ( ), "r" ) );

            for ( int dim = 0; dim < tdim ( ) + 1; ++dim )
            {
                std::vector<std::string> names ( 0 );
                std::stringstream attr_names;
                attr_names << "attribute_names";
                if ( dim < 10 )
                    attr_names << "0";
                attr_names << dim;
                read_array_parallel ( group_ptr, attr_names.str ( ), names, comm );
                if ( names[0] == "__empty__" ) continue;
                for ( std::vector<std::string>::const_iterator it = names.begin ( ); it != names.end ( ); ++it )
                {
                    char att_type = ( *it )[0];
                    std::string attribute_name = ( *it ).substr ( 1 );
                    AttributePtr attribute;
                    if ( att_type == 'i' )
                    {
                        std::vector<int> int_att ( 0 );
                        read_array_parallel<int>( group_ptr, attribute_name, int_att, comm );
                        attribute = AttributePtr ( new IntAttribute ( int_att ) );
                    }
                    else if ( att_type == 'd' )
                    {
                        std::vector<double> double_att ( 0 );
                        read_array_parallel<double>( group_ptr, attribute_name, double_att, comm );
                        attribute = AttributePtr ( new DoubleAttribute ( double_att ) );
                    }
                    else if ( att_type == 'b' )
                    {
                        std::vector<int> temp_bool_att ( 0 );
                        read_array_parallel<int>( group_ptr, attribute_name, temp_bool_att, comm );
                        std::vector<bool> bool_att ( temp_bool_att.size ( ) );
                        for ( int i = 0; i < temp_bool_att.size ( ); ++i )
                        {
                            bool_att[i] = ( temp_bool_att[i] == 1 );
                        }
                        attribute = AttributePtr ( new BoolAttribute ( bool_att ) );
                    }
                    add_attribute ( attribute_name, dim, attribute );
                }
            }

            // Load period
            std::vector<double> conv_period;
            read_array_parallel<double>( group_ptr, "period", conv_period, comm );
            if ( conv_period.size ( ) == 1 ) //case period is empty
                conv_period.pop_back ( );
            assert ( conv_period.size ( ) % 4 == 0 );
            int period_size = conv_period.size ( ) / 4;
            period_.resize ( period_size );
            for ( int k = 0; k < period_size; ++k )
            {
                double master = conv_period[4 * k + 0];
                double slave = conv_period[4 * k + 1];
                double h = conv_period[4 * k + 2];
                int index = ( int ) conv_period[4 * k + 3];
                period_[k] = MasterSlave ( master, slave, h, index );
            }
#    else
            throw "HiFlow was not compiled with HDF5 support!\n";
#    endif
        }
    } // namespace mesh
} // namespace hiflow

#endif
